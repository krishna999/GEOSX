/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "VTKPolyDataWriterInterface.hpp"
#include "dataRepository/Group.hpp"

// TPL includes
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkExtentTranslator.h>

// System includes
#include <unordered_set>
#include <sys/stat.h>



namespace geosx
{
namespace vtk
{


/*!
 * @brief Get the DataSet file path
 * @details the DataSet file path is the path to the .vtu mesh file
 * @param[in] er The ElementRegion
 * @param[in] time the time-step
 * @param[in] rank the rank to be written
 */
string GetDataSetFilePath( ElementRegionBaseABC const & er, double time, int rank )
{
  return std::to_string( time ) + "/" + stringutilities::PadValue( rank, std::to_string( MpiWrapper::commSize() ).size() ) + "_" + er.getNameMock() + ".vtu";
}

const std::map< string, int > geosx2VTKCellTypes =
{
  { "C3D4", VTK_TETRA },
  { "C3D8", VTK_HEXAHEDRON },
  { "C3D6", VTK_WEDGE },
  { "C3D5", VTK_PYRAMID }
};

/*!
 * @brief Gets the VTK cell identifier
 * @param[in] elementType the type of the element (using the abaqus nomenclature)
 * @return the VTK cell identifier
 */
int ToVTKCellType( const string & elementType )
{
  int vtkIdentifier = VTK_EMPTY_CELL;
  try
  {
    vtkIdentifier = geosx2VTKCellTypes.at( elementType );
  }
  catch( const std::out_of_range & outOfRange )
  {
    GEOSX_ERROR( "Element type " << elementType << " not recognized for VTK output " );
  }
  return vtkIdentifier;
}

void VTKPolyDataWriterInterface::gatherNbElementsInRegion( ElementRegionBaseABC const & er,
                                                           array1d< localIndex > & nbElemsInRegion ) const
{
  localIndex nbElems = er.nElementsMock();
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  nbElemsInRegion.resize( mpiSize );
  MpiWrapper::gather( &nbElems, 1, nbElemsInRegion.data(), 1, 0, MPI_COMM_GEOSX );
}

VTKPolyDataWriterInterface::VTKPolyDataWriterInterface( string const & outputName ):
  m_outputFolder( outputName ),
  m_pvd( outputName + ".pvd" ),
  m_previousCycle( -1 )
{
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    mode_t mode = 0733;
    int errorCode = mkdir( outputName.c_str(), mode );
    if( errorCode == -1 )
    {
      if( errno == EEXIST )
      {
        GEOSX_WARNING( "Path " << outputName << " already exists from a previous simulation" );
      }
      else
      {
        GEOSX_ERROR( "Fail to create main directory " << outputName << " for the VTK output" );
      }
    }
  }
  MpiWrapper::barrier();
}

vtkSmartPointer< vtkPoints >  VTKPolyDataWriterInterface::getVtkPoints( NodeManagerABC const & nodeManager ) const
{
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  points->SetNumberOfPoints( nodeManager.nPoints() );
  points->SetNumberOfPoints( nodeManager.nPoints() );
  auto connectivity = nodeManager.referencePosition();
  for( localIndex v = 0; v < nodeManager.nPoints(); v++ )
  {
    points->SetPoint( v, connectivity[v][0], connectivity[v][1], connectivity[v][2] );
  }
  return points;
}

std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >VTKPolyDataWriterInterface::getWell( WellElementSubRegionABC const & esr,
                                                                                                               NodeManagerABC const & nodeManager ) const
{
  // some notes about WellElementSubRegion:
  // - if the well represented by this subRegion is not on this rank, esr.size() = 0
  // - otherwise, esr.size() is equal to the number of well elements of the well on this rank
  // Each well element has two nodes, shared with the previous and next well elements, respectively
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  // if esr.size() == 0, we set the number of points and cells to zero
  // if not, we set the number of points to esr.size()+1 and the number of cells to esr.size()
  localIndex const numPoints = esr.nCells() > 0 ? esr.nCells() + 1 : 0;
  points->SetNumberOfPoints( numPoints );
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( esr.nCells() );
  localIndex const numberOfNodesPerElement = esr.nNodesPerElementMock();
  GEOSX_ERROR_IF_NE( numberOfNodesPerElement, 2 );
  std::vector< vtkIdType > connectivity( numberOfNodesPerElement );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const referencePosition = nodeManager.referencePosition();

  // note that if esr.size() == 0, we don't have any point or cell to add below and we just return

  for( localIndex edge = 0; edge < esr.nCells(); edge++ )
  {
    localIndex const firstPoint = esr.nodeList()[edge][0];
    auto point = referencePosition[firstPoint];
    points->SetPoint( edge, point[0], point[1], point[2] );
    connectivity[0] = edge;
    connectivity[1] = edge+1; // We can do that because of the pattern in which the wells are stored
    cellsArray->InsertNextCell( numberOfNodesPerElement, connectivity.data() );
  }

  if( esr.nCells() > 0 )
  {
    localIndex const lastPoint = esr.nodeList()[ esr.nCells() -1  ][1];
    auto point = referencePosition[lastPoint];
    points->SetPoint( esr.nCells(), point[0], point[1], point[2] );
  }

  return std::make_pair( points, cellsArray );
}

std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
VTKPolyDataWriterInterface::getSurface( FaceElementSubRegionABC const & esr,
                                        NodeManagerABC const & nodeManager ) const
{
  // Get unique node set composing the surface
  auto & nodeListPerElement = esr.nodeList();
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( esr.nCells() );
  std::unordered_map< localIndex, localIndex > geosx2VTKIndexing;
  geosx2VTKIndexing.reserve( esr.nCells() * esr.nNodesPerElementMock() );
  localIndex nodeIndexInVTK = 0;
  std::vector< vtkIdType > connectivity( esr.nNodesPerElementMock() );
  std::vector< int >  vtkOrdering = esr.getVTKNodeOrderingMock();

  for( localIndex ei = 0; ei < esr.nCells(); ei++ )
  {
    auto const & elem = nodeListPerElement[ei];
    for( localIndex i = 0; i < elem.size(); i++ )
    {
      auto const & VTKIndexPos = geosx2VTKIndexing.find( elem[vtkOrdering[i]] );
      if( VTKIndexPos == geosx2VTKIndexing.end() )
      {
        connectivity[i] = geosx2VTKIndexing[elem[vtkOrdering[i]]] = nodeIndexInVTK++;
      }
      else
      {
        connectivity[i] = VTKIndexPos->second;
      }
    }
    cellsArray->InsertNextCell( elem.size(), connectivity.data() );
  }

  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  points->SetNumberOfPoints( geosx2VTKIndexing.size() );
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const referencePosition = nodeManager.referencePosition();

  for( auto nodeIndex: geosx2VTKIndexing )
  {
    auto point = referencePosition[nodeIndex.first];
    points->SetPoint( nodeIndex.second, point[0], point[1], point[2] );
  }

  return std::make_pair( points, cellsArray );
}
std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
VTKPolyDataWriterInterface::getEmbeddedSurface( EmbeddedSurfaceSubRegionABC const & esr,
                                                NodeManagerABC const & nodeManager ) const
{
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & intersectionPoints = nodeManager.embSurfNodesPosition();

  points->SetNumberOfPoints( intersectionPoints.size( 0 ) );
  for( localIndex pointIndex = 0; pointIndex < intersectionPoints.size( 0 ); pointIndex++ )
  {
    points->SetPoint( pointIndex, intersectionPoints[pointIndex][0], intersectionPoints[pointIndex][1], intersectionPoints[pointIndex][2] );
  }

  EmbeddedSurfaceSubRegion::NodeMapType const & toNodesMap = esr.nodeList();
  std::vector< vtkIdType > connectivity( 10 );
  for( localIndex cellIndex = 0; cellIndex < esr.nCells(); cellIndex++ )
  {
    connectivity.resize( toNodesMap.sizeOfArray( cellIndex ) );
    for( std::size_t i = 0; i < connectivity.size(); ++i )
    {
      connectivity[i] = esr.nodeList( cellIndex, i );
    }
    cellsArray->InsertNextCell( connectivity.size(), connectivity.data() );
  }

  return std::make_pair( points, cellsArray );
}

std::pair< std::vector< int >, vtkSmartPointer< vtkCellArray > >
VTKPolyDataWriterInterface::getVtkCells( CellElementRegionABC const & er ) const
{
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( er.getNElementsInCellElementSubRegion() );
  std::vector< int > cellType;
  cellType.reserve( er.getNElementsInCellElementSubRegion() );
  for( const CellElementSubRegionABC & esr: er.getCellElementSubRegions() )
  {
    std::vector< vtkIdType > connectivity( esr.nNodesPerElementMock() );
    std::vector< int > vtkOrdering = esr.getVTKNodeOrderingMock();
    int vtkCellType = ToVTKCellType( esr.getElementTypeString() );
    for( localIndex c = 0; c < esr.size(); c++ ) // FIXME size of what?
    {
      for( std::size_t i = 0; i < connectivity.size(); i++ )
      {
        connectivity[i] = esr.nodeList( c, vtkOrdering[i] );
      }
      cellType.push_back( vtkCellType );
      cellsArray->InsertNextCell( esr.nNodesPerElementMock(), connectivity.data() );
    }
  }
  return std::make_pair( cellType, cellsArray );
}

void VTKPolyDataWriterInterface::writeField( WrapperBase const & wrapperBase,
                                             vtkSmartPointer< VTKGEOSXData > data,
                                             localIndex size, localIndex & count ) const
{
  std::type_info const & typeID = wrapperBase.getTypeId();
  if( typeID==typeid(r1_array) )
  {
    // We need a special case for the R1 array
    // Because those array are not stored the same way than the classical array2d which is
    // preferable to store a vector on each element.
    data->SetNumberOfComponents( 3 );
  }
  rtTypes::applyArrayTypeLambda2( rtTypes::typeID( typeID ),
                                  true,
                                  [&]( auto array, auto GEOSX_UNUSED_PARAM( Type ) )->void
  {
    typedef decltype( array ) arrayType;
    Wrapper< arrayType > const & wrapperT = Wrapper< arrayType >::cast( wrapperBase );
    traits::ViewTypeConst< arrayType > const sourceArray = wrapperT.reference().toViewConst();
    if( typeID!=typeid(r1_array) )
    {
      integer nbOfComponents = 1;
      for( localIndex i = 1; i < arrayType::NDIM; i++ )
      {
        nbOfComponents = nbOfComponents * sourceArray.size( i );
      }
      data->SetNumberOfComponents( nbOfComponents );
    }
    for( localIndex i = 0; i < size; i++ )
    {
      LvArray::forValuesInSlice( sourceArray[i], [&]( auto const & value )
      {
        data->customInsertValue( count++, value );
      } );
    }
  } );

}

void VTKPolyDataWriterInterface::writeNodeFields( vtkSmartPointer< vtkPointData > const pointdata,
                                                  NodeManagerABC const & nodeManager ) const
{
  for( auto const & wrapper : nodeManager.getWrappers() )
  {
//    auto const & wrapper = *wrapperIter.second;
    if( wrapper.getPlotLevel() <= m_plotLevel )
    {
      vtkSmartPointer< VTKGEOSXData > data = VTKGEOSXData::New();
      data->SetNumberOfValues( nodeManager.nPoints() );
      data->SetName( wrapper.getName().c_str() );;
      localIndex count = 0;
      writeField( wrapper, data, nodeManager.nPoints(), count );
      pointdata->AddArray( data );
    }
  }
}

//template< class SUBREGION >
//void VTKPolyDataWriterInterface::writeElementFields( vtkSmartPointer< vtkCellData > const celldata,
//                                                     ElementRegionBaseABC const & er ) const
//{
//  std::unordered_set< string > allFields;
//  er.forElementSubRegions< SUBREGION >( [&]( auto const & esr )
//  {
//    for( auto const & wrapperIter : esr.wrappers() )
//    {
//      auto const * const wrapper = wrapperIter.second;
//      if( wrapper->getPlotLevel() <= m_plotLevel )
//      {
//        allFields.insert( wrapperIter.first );
//      }
//    }
//  } );
//
//  for( auto const & field : allFields )
//  {
//    vtkSmartPointer< VTKGEOSXData > data = VTKGEOSXData::New();
//    data->SetNumberOfValues( er.getNumberOfElements< SUBREGION >() );
//    data->SetName( field.c_str() );
//
//    localIndex count = 0;
//    er.forElementSubRegions< SUBREGION >( [&]( auto const & esr )
//    {
//      auto const & wrapper = *esr.getWrapperBase( field );
//      WriteField( wrapper, data, esr.size(), count );
//    } );
//    celldata->AddArray( data );
//  }
//}

//template< class ELEMENT_REGION >
//void VTKPolyDataWriterInterface::MyWriteElementFields( vtkSmartPointer< vtkCellData > const celldata,
//                                                       ELEMENT_REGION const & er ) const
//{
//  std::unordered_set< string > allFields;
//  for( const auto & esr: er.getElementSubRegionsSpec() )
//  {
//    for( auto const & wrapperIter : esr.wrappers() )
//    {
//      auto const * const wrapper = wrapperIter.second;
//      if( wrapper->getPlotLevel() <= m_plotLevel )
//      {
//        allFields.insert( wrapperIter.first );
//      }
//    }
//  }
//
//  for( auto const & field : allFields )
//  {
//    vtkSmartPointer< VTKGEOSXData > data = VTKGEOSXData::New();
//    data->SetNumberOfValues( er.getSpecElementNbr() );
//    data->SetName( field.c_str() );
//
//    localIndex count = 0;
//    for( const auto & esr: er.getElementSubRegionsSpec() )
//    {
//      auto const & wrapper = *esr.getWrapperBase( field );
//      WriteField( wrapper, data, esr.size(), count );
//    }
//
//    celldata->AddArray( data );
//  }
//}

template< typename CONTAINER_TYPE, typename CAST_TYPE >
static bool MyMatch( CONTAINER_TYPE & container )
{
  // FIXME facotrize as private static member function...
  using T = std::conditional_t< std::is_const< CONTAINER_TYPE >::value, CAST_TYPE const, CAST_TYPE >;
  T * const castedContainer = dynamic_cast< T * >( &container );

  return castedContainer != nullptr;
//  if( castedContainer != nullptr )
//  {
//    return true;
//  }
//
//  return false;
}

//template< typename CONTAINER_TYPE, typename T0, typename T1, typename ... CAST_TYPES >
//struct MyPredicate
//{
//public:
template< typename CONTAINER_TYPE, typename T0, typename T1, typename ... CAST_TYPES >
bool MyMatch( const CONTAINER_TYPE & container )
{
  using T = std::conditional_t< std::is_const< CONTAINER_TYPE >::value, T0 const, T0 >;
  T * const castedContainer = dynamic_cast< T * >( &container );

  if( castedContainer != nullptr )
  {
    return true;
  }

  return MyMatch< T1, CAST_TYPES... >( container );
}

template< typename CAST_TYPE >
struct TestMatch
{
  template< typename CONTAINER_TYPE >
  static bool test(CONTAINER_TYPE & container)
  {
    using T = std::conditional_t< std::is_const< CONTAINER_TYPE >::value, CAST_TYPE const, CAST_TYPE >;
    T * const castedContainer = dynamic_cast< T * >( &container );

    return castedContainer != nullptr;
  }
};

template< class PREDICATE >
void VTKPolyDataWriterInterface::writeElementFieldsPredicate( vtkSmartPointer< vtkCellData > const celldata, ElementRegionBaseABC const & er, PREDICATE p ) const
{
  std::unordered_set< string > allFields;
  for( const auto & esr: er.getElementSubRegions() )
  {
    if( not p( esr ) ) continue;

    for( auto const & wrapperIter : esr.wrappers() )
    {
      auto const * const wrapper = wrapperIter.second;
      if( wrapper->getPlotLevel() <= m_plotLevel )
      {
        allFields.insert( wrapperIter.first );
      }
    }
  }

  for( auto const & field : allFields )
  {
    vtkSmartPointer< VTKGEOSXData > data = VTKGEOSXData::New();
    data->SetNumberOfValues( er.getNumberOfElements< SUBREGION >() );
    data->SetName( field.c_str() );

    localIndex count = 0;
    for( const auto & esr: er.getElementSubRegions() ) //  FIXME use std::ranges::views::filter here
    {
      if( not p( esr ) ) continue;

      auto const & wrapper = *esr.getWrapperBase( field );
      writeField( wrapper, data, esr.get().nCells(), count );
    }
    celldata->AddArray( data );
  }
}

void VTKPolyDataWriterInterface::writeCellElementRegions( real64 time,
                                                          ElementRegionManagerABC const & elemManager,
                                                          NodeManagerABC const & nodeManager ) const
{
//  auto pred = [](const ElementSubRegionBaseABC & esr)
//  {
//    return MyMatch< ElementSubRegionBaseABC, CellElementSubRegionABC >( esr );
//  };
//  auto predicate = MyMatch< const ElementSubRegionBaseABC, CellElementSubRegionABC >; // FIXME GOOD?
  auto predicate = []( const ElementSubRegionBaseABC & esr ) -> bool // FIXME this is mandatory to infer the PREDICATE type (because we need to infer the CONTAINER_TYPE)
  {
    return TestMatch< CellElementSubRegionABC >::test( esr );
  };
  for (const CellElementRegionABC & er: elemManager.getCellElementRegions() )
  {
    if( er.getNElementsInCellElementSubRegion() != 0 ) // FIXME was it relevant to have template args other than CellElementSubRegion?`
    {
      vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
      auto VTKPoints = getVtkPoints( nodeManager );
      ug->SetPoints( VTKPoints );
      auto VTKCells = getVtkCells( er );
      ug->SetCells( VTKCells.first.data(), VTKCells.second );
//      writeElementFields< CellElementSubRegionABC >( ug->GetCellData(), er );
      writeElementFieldsPredicate( ug->GetCellData(), er, predicate );
//      writeElementFieldsPredicate( ug->GetCellData(), er, MyMatch< const ElementSubRegionBaseABC, CellElementSubRegionABC > );
      writeNodeFields( ug->GetPointData(), nodeManager );
      writeUnstructuredGrid( ug, time, er.getNameMock() );
    }
  }
}

void VTKPolyDataWriterInterface::writeWellElementRegions( real64 time, ElementRegionManagerABC const & elemManager,
                                                          NodeManagerABC const & nodeManager ) const
{
  auto predicate = []( const ElementSubRegionBaseABC & esr ) -> bool
  {
    return TestMatch< WellElementSubRegionABC >::test( esr );
  };

  for (const WellElementRegionABC & er: elemManager.getWellElementRegions() )
  {
    auto esr = er.getWellElementSubRegion();
    vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
    auto VTKWell = getWell( *esr, nodeManager );
    ug->SetPoints( VTKWell.first );
    ug->SetCells( VTK_LINE, VTKWell.second );
//    writeWellElementRegions< WellElementSubRegionABC >( ug->GetCellData(), er );
//    WriteElementFieldsPredicate( ug->GetCellData(), er, MyMatch< const ElementSubRegionBaseABC, WellElementSubRegionABC > );
    writeElementFields< WellElementSubRegion >( ug->GetCellData(), er );
    writeUnstructuredGrid( ug, time, er.getNameMock() );
  } );
}

void VTKPolyDataWriterInterface::writeSurfaceElementRegions( real64 time,
                                                             ElementRegionManagerABC const & elemManager,
                                                             NodeManagerABC const & nodeManager ) const
{
  for (const SurfaceElementRegionABC & er: elemManager.getSurfaceElementRegions() )
  {
    vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
    if( er.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::embeddedElement )
    {
      auto predicate = []( const ElementSubRegionBaseABC & esr ) -> bool
      {
        return TestMatch< EmbeddedSurfaceSubRegionABC >::test( esr );
      };

      auto esr = er.getEmbeddedSurfaceSubRegion();

      auto VTKSurface = getEmbeddedSurface( *esr, nodeManager );
      ug->SetPoints( VTKSurface.first );
      ug->SetCells( VTK_POLYGON, VTKSurface.second );

//      writeElementFields< EmbeddedSurfaceSubRegionABC >( ug->GetCellData(), er );
      writeElementFieldsPredicate( ug->GetCellData(), er, predicate );
    }
    else if( er.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::faceElement )
    {
      auto predicate = []( const ElementSubRegionBaseABC & esr ) -> bool
      {
        return TestMatch< FaceElementSubRegionABC >::test( esr );
      };

      auto esr = er.getFaceElementSubRegion();

      auto VTKSurface = getSurface( *esr, nodeManager );

      ug->SetPoints( VTKSurface.first );
      if( esr->nNodesPerElementMock() == 8 )
      {
        ug->SetCells( VTK_HEXAHEDRON, VTKSurface.second );
      }
      else if( esr->nNodesPerElementMock() == 6 )
      {
        ug->SetCells( VTK_WEDGE, VTKSurface.second );
      }
      else
      {
        GEOSX_ERROR( "Elements with " << esr->nNodesPerElementMock() << " nodes can't be output "
                                      << "in the FaceElementRegion " << er.getNameMock() );
      }
//      writeElementFields< FaceElementSubRegionABC >( ug->GetCellData(), er );
      writeElementFieldsPredicate( ug->GetCellData(), er, predicate );
    }
    writeUnstructuredGrid( ug, time, er.getNameMock() );
  }
}

void VTKPolyDataWriterInterface::writeVtmFile( real64 time,
                                               ElementRegionManagerABC const & elemManager,
                                               VTKVTMWriter const & vtmWriter ) const
{
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  auto writeSubBlocks = [&]( ElementRegionBaseABC const & er )
  {
    array1d< localIndex > nbElemsInRegion;
    gatherNbElementsInRegion( er, nbElemsInRegion );
    vtmWriter.addSubBlock( er.getCatalogNameMock(), er.getNameMock() );
    for( int i = 0; i < mpiSize; i++ )
    {
      string const dataSetFile = GetDataSetFilePath( er, time, i );
      if( nbElemsInRegion[i] > 0 )
      {
        if( mpiRank == 0 )
        {
          vtmWriter.addDataToSubBlock( er.getCatalogNameMock(), er.getNameMock(), dataSetFile, i );
        }
      }
    }
  };
  if( mpiRank == 0 )
  {
    vtmWriter.addBlock( CellElementRegion::catalogName() );
    vtmWriter.addBlock( WellElementRegion::catalogName() );
    vtmWriter.addBlock( SurfaceElementRegion::catalogName() );
  }

  for( const CellElementRegionABC & cer: elemManager.getCellElementRegions() )
    writeSubBlocks( cer );
  for( const WellElementRegionABC & wer: elemManager.getWellElementRegions() )
    writeSubBlocks( wer );
  for( const SurfaceElementRegionABC & ser: elemManager.getSurfaceElementRegions() )
    writeSubBlocks( ser );

  if( mpiRank == 0 )
  {
    vtmWriter.save();
  }
}

void VTKPolyDataWriterInterface::createTimeStepSubFolder( real64 time ) const
{
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    mode_t mode = 0773;
    string const timeStepSubFolder = getTimeStepSubFolder( time );
    int errorCode = mkdir( timeStepSubFolder.c_str(), mode );
    if( errorCode == -1 )
    {
      if( errno == EEXIST )
      {
        GEOSX_WARNING( "Path " << timeStepSubFolder << " already exists from a previous simulation" );
      }
      else
      {
        GEOSX_ERROR( "Fail to create the timestep directory " << timeStepSubFolder << " for the VTK output" );
      }
    }
  }
  MpiWrapper::barrier();
}

void VTKPolyDataWriterInterface::writeUnstructuredGrid( vtkSmartPointer< vtkUnstructuredGrid > ug,
                                                        double time,
                                                        string const & name ) const
{
  string timeStepSubFolder = VTKPolyDataWriterInterface::getTimeStepSubFolder( time );
  vtkSmartPointer< vtkXMLUnstructuredGridWriter > vtuWriter =vtkXMLUnstructuredGridWriter::New();
  vtuWriter->SetInputData( ug );
  string vtuFilePath = timeStepSubFolder + "/" +
                       stringutilities::PadValue( MpiWrapper::commRank(), std::to_string( MpiWrapper::commSize() ).size() ) +"_" + name + ".vtu";
  vtuWriter->SetFileName( vtuFilePath.c_str() );
  if( m_outputMode == VTKOutputMode::BINARY )
  {
    vtuWriter->SetDataModeToBinary();
  }
  else if( m_outputMode == VTKOutputMode::ASCII )
  {
    vtuWriter->SetDataModeToAscii();
  }
  vtuWriter->Write();
}

string VTKPolyDataWriterInterface::getTimeStepSubFolder( real64 time ) const
{
  return m_outputFolder + "/" + std::to_string( time );
}

void VTKPolyDataWriterInterface::write( real64 time, integer cycle, DomainPartition const & domain )
{
  createTimeStepSubFolder( time );
  ElementRegionManager const & elemManager = *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  NodeManager const & nodeManager = *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getNodeManager();
  writeCellElementRegions( time, elemManager, nodeManager );
  writeWellElementRegions( time, elemManager, nodeManager );
  writeSurfaceElementRegions( time, elemManager, nodeManager );
  string vtmPath = getTimeStepSubFolder( time ) + ".vtm";
  VTKVTMWriter vtmWriter( vtmPath );
  writeVtmFile( time, elemManager, vtmWriter );
  if( cycle != m_previousCycle )
  {
    m_pvd.addData( time, vtmPath );
    m_pvd.save();
  }
  m_previousCycle = cycle;
}

void VTKPolyDataWriterInterface::write( real64 time,
                                        integer cycle,
                                        ElementRegionManagerABC const & elemManager,
                                        NodeManagerABC const & nodeManager )
{
  createTimeStepSubFolder( time );
  writeCellElementRegions( time, elemManager, nodeManager );
  writeWellElementRegions( time, elemManager, nodeManager );
  writeSurfaceElementRegions( time, elemManager, nodeManager );
  string vtmPath = getTimeStepSubFolder( time ) + ".vtm";
  VTKVTMWriter vtmWriter( vtmPath );
  writeVtmFile( time, elemManager, vtmWriter );
  if( cycle != m_previousCycle )
  {
    m_pvd.addData( time, vtmPath );
    m_pvd.save();
  }
  m_previousCycle = cycle;
}

}
}
