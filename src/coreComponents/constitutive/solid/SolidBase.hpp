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


/**
 * @file SolidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
namespace constitutive
{

/**
 * @brief Base class for all solid constitutive kernel wrapper classes.
 *
 * The responsibility of this base is to:
 *
 * 1) Contain views to state and parameter data for solid models.
 * 2) Specify an interface for state update functions.
 *
 * In general, the ArrayView data in the wrapper is specified to be of type
 * "arrayView<T> const" or "arrayView<T const> const". The "const-ness"
 * of the data indicates whether it is a parameter" or a state variable,
 * with the parameters being "T const" and state variables being "T".
 *
 * @note
 * If an allocation occurs on the underlying Array after a KernelWrapper is created,
 * then the ArrayView members of that KernelWrapper are silently invalid.
 */
class SolidBaseUpdates
{
protected:
  /**
   * @brief constructor
   * @param[in] newStress The new stress data from the constitutive model class.
   * @param[in] oldStress The old stress data from the constitutive model class.
   */
  SolidBaseUpdates( arrayView3d< real64, solid::STRESS_USD > const & newStress,
                    arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    m_newStress( newStress ),
    m_oldStress( oldStress )
  {}

  /// Deleted default constructor
  SolidBaseUpdates() = delete;

  /**
   * @brief Copy Constructor
   * @param source Object to copy
   */
  SolidBaseUpdates( SolidBaseUpdates const & source ) = default;

  /**
   * @brief Move Constructor
   * @param source Object to move resources from
   */
  SolidBaseUpdates( SolidBaseUpdates && source ) = default;

  /// Deleted copy assignment operator
  SolidBaseUpdates & operator=( SolidBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  SolidBaseUpdates & operator=( SolidBaseUpdates && ) =  delete;

  /**
   * @brief Helper to save point stress back to m_newStress array
   *
   * This is mostly defined for improving code readability.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] stress Stress to be save to m_newStress[k][q]
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void saveStress( localIndex const k,
                   localIndex const q,
                   real64 const ( &stress )[6] ) const
  {
    LvArray::tensorOps::copy< 6 >( m_newStress[k][q], stress );
  }


public:

  /// A reference the current material stress at quadrature points.
  arrayView3d< real64, solid::STRESS_USD > const m_newStress;

  /// A reference the previous material stress at quadrature points.
  arrayView3d< real64, solid::STRESS_USD > const m_oldStress;

  /**
   * @name Update Interfaces: Stress and Stiffness
   *
   * We define a variety of interfaces for constitutive models using different
   * strain theories.  Derived classes only need to implement the subset of interfaces
   * most relevant to them.
   *
   * This group of interfaces returns stress and stiffness simultaneously, and
   * are most useful for implicit finite element formulations.
   */
  ///@{

  /**
   * @brief Small strain update.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] strainIncrement Strain increment in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New tangent stiffness value
   */
  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  real64 ( & stiffness )[6][6] ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( strainIncrement );
    GEOSX_UNUSED_VAR( stress );
    GEOSX_UNUSED_VAR( stiffness );
    GEOSX_ERROR( "smallStrainUpdate() not implemented for this model" );
  }

  /**
   * @brief Small strain, stateless update.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] totalStrain Total strain in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New tangent stiffness value
   */
  GEOSX_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( & stress )[6],
                                         real64 ( & stiffness )[6][6] ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( totalStrain );
    GEOSX_UNUSED_VAR( stress );
    GEOSX_UNUSED_VAR( stiffness );
    GEOSX_ERROR( "smallStrainNoStateUpdate() not implemented for this model" );
  }

  /**
   * @brief Hypo update (small strain, large rotation).
   *
   * The base class uses a call to the small strain update, followed by
   * a rotation correction using the Hughes-Winget incrementally objective
   * algorithm.  One can imagine the material deforming in small strain, but in
   * a reference frame that rotates to track the body's local rotation.  This
   * provides a convenient way to extend small strain models to a finite rotation
   * regime.  From the assumption of small deformations, the Cauchy stress and
   * Kirchoff stress are approximately equal (det F ~ 1).
   *
   * Note that if the derived class has tensorial state variables (beyond the
   * stress itself) care must be taken to rotate these as well.
   *
   * We use a post-rotation of the new stress (as opposed to pre-rotation of the old stress)
   * as we don't have to unrotate the old stress in the event a rewind is required.
   * We do not post-rotate the stiffness tensor, which is an approximation, but
   * should be sufficient for small rotation increments.
   *
   * This method does not work for anisotropic properties or yield functions
   * (without some care) and a co-rotational formulation should be considered instead.
   *
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] Ddt The incremental deformation tensor (rate of deformation tensor * dt)
   * @param[in] Rot The incremental rotation tensor
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New stiffness value
   */
  GEOSX_HOST_DEVICE
  virtual void hypoUpdate( localIndex const k,
                           localIndex const q,
                           real64 const ( &Ddt )[6],
                           real64 const ( &Rot )[3][3],
                           real64 ( & stress )[6],
                           real64 ( & stiffness )[6][6] ) const
  {
    smallStrainUpdate( k, q, Ddt, stress, stiffness );

    real64 temp[6] = { 0 };
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( temp, Rot, m_newStress[ k ][ q ] );
    LvArray::tensorOps::copy< 6 >( stress, temp );
    saveStress( k, q, stress );
  }

  /**
   * @brief Hyper update (large deformation).
   *
   * This version provides an interface for fully-general, large-deformation
   * models.  The input strain measure is the deformation gradient minus the identity.
   * The output is the Cauchy stress (true stress in the deformed configuration) to be
   * consistent with the previous interfaces.  The stiffness is similarly the
   * stiffness in the deformed configuration.
   *
   * @note: This interface is currently a placeholder and has not been extensively
   * tested.
   *
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] FmI Deformation gradient minus identity (F-I)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New stiffness value
   */
  // TODO: confirm stress and strain measures we want to use
  GEOSX_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( &FminusI )[3][3],
                            real64 ( & stress )[6],
                            real64 ( & stiffness )[6][6] ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( FminusI );
    GEOSX_UNUSED_VAR( stress );
    GEOSX_UNUSED_VAR( stiffness );
    GEOSX_ERROR( "hyperUpdate() not implemented for this model" );
  }

  ///@}
  /**
   * @name Update Interfaces: Stress-Only
   *
   * We define a variety of interfaces for constitutive models using different
   * strain theories.  Derived classes only need to implement the subset of interfaces
   * most relevant to them.
   *
   * This group of interfaces returns stress only, with no stiffness, and
   * are most useful for explicit finite element formulations.
   *
   * @note
   * The base class versions implement a naive call to the versions
   * of the updates that returns both stress and stiffness, and then discards
   * the stiffness.  Derived classes can implement an optimized version to
   * avoid extranenous work, but we delegate this detail to them.
   */
  ///@{

  /**
   * @brief Small strain update, returning only stress.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] strainIncrement Strain increment in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New tangent stiffness value
   */
  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( & stress )[6] ) const
  {
    real64 discard[6][6];
    smallStrainUpdate( k, q, strainIncrement, stress, discard );
  }


  /**
   * @brief Small strain, stateless update, returning only stress.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] totalStrain total strain in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New tangent stiffness value
   */
  GEOSX_HOST_DEVICE
  virtual void smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( & stress )[6] ) const
  {
    real64 discard[6][6];
    smallStrainNoStateUpdate( k, q, totalStrain, stress, discard );
  }

  /**
   * @brief Hypo update, returning only stress
   *
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] Ddt The incremental deformation tensor (rate of deformation tensor * dt)
   * @param[in] Rot The incremental rotation tensor
   * @param[out] stress New stress value (Cauchy stress)
   */
  GEOSX_HOST_DEVICE
  virtual void hypoUpdate_StressOnly( localIndex const k,
                                      localIndex const q,
                                      real64 const ( &Ddt )[6],
                                      real64 const ( &Rot )[3][3],
                                      real64 ( & stress )[6] ) const
  {
    real64 discard[6][6];
    hypoUpdate( k, q, Ddt, Rot, stress, discard );
  }

  /**
   * @brief Hyper update, returning only stresses.
   *
   * @param[in] k The element index.
   * @param[in] q The quadrature point index.
   * @param[in] FmI Deformation gradient minus identity (F-I)
   * @param[out] stress New stress value (Cauchy stress)
   */
  GEOSX_HOST_DEVICE
  virtual void hyperUpdate_StressOnly( localIndex const k,
                                       localIndex const q,
                                       real64 const ( &FminusI )[3][3],
                                       real64 ( & stress )[6] ) const
  {
    real64 discard[6][6];
    hyperUpdate( k, q, FminusI, stress, discard );
  }

  ///@}


  /**
   * @brief Return the current elastic strain at a given material point (small-strain interface)
   *
   * @param k the element inex
   * @param q the quadrature index
   * @param elasticStrain Current elastic strain
   */
  GEOSX_HOST_DEVICE
  virtual void getElasticStrain( localIndex const k,
                                 localIndex const q,
                                 real64 ( & elasticStrain )[6] ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( elasticStrain );
    GEOSX_ERROR( "getElasticStrain() not implemented for this model" );
  }

  /**
   * @brief Return the strain energy density at a given material point
   *
   * @param k the element inex
   * @param q the quadrature index
   * @return Strain energy density
   */
  GEOSX_HOST_DEVICE
  virtual real64 getStrainEnergyDensity( localIndex const k,
                                         localIndex const q ) const
  {
    auto const & stress = m_newStress[k][q];

    real64 strain[6];
    getElasticStrain( k, q, strain );

    real64 energy = 0;

    for( localIndex i=0; i<6; ++i )
    {
      energy += stress[i]*strain[i];  // contraction sigma:epsilon
    }
    energy *= 0.5;

    GEOSX_ASSERT_MSG( energy >= 0.0, "negative strain energy density detected" );

    return energy;
  }

  /**
   * @brief Return the stiffness at a given element (small-strain interface)
   *
   * @note If the material model has a strain-dependent material stiffness (e.g.
   * any plasticity, damage, or nonlinear elastic model) then this interface will
   * not work.  Users should instead use one of the interfaces where a strain
   * tensor is provided as input.
   *
   * @note Given the limitations above, this function may be removed from the
   * public interface in the future.  Direct use in physics
   * solvers is discouraged.
   *
   * @param k the element number
   * @param stiffness the stiffness array
   */
  GEOSX_HOST_DEVICE
  virtual void getElasticStiffness( localIndex const k, real64 ( & stiffness )[6][6] ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( stiffness );
    GEOSX_ERROR( "getElasticStiffness() not implemented for this model" );
  }

};


/**
 * @class SolidBase
 * This class serves as the base class for solid constitutive models.
 */
class SolidBase : public constitutive::ConstitutiveBase
{
public:
  /**
   * @brief Constructor
   * @param name Name of the SolidBase object in the repository.
   * @param parent The parent group of the SolidBase object.
   */
  SolidBase( string const & name,
             Group * const parent );

  /**
   * Destructor
   */
  virtual ~SolidBase() override;

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto stressString = "stress";
    static constexpr auto oldStressString = "oldStress";
    static constexpr auto densityString  = "density";
    static constexpr auto defaultDensityString  = "defaultDensity";
  };


  /// Number of elements storing solid data
  localIndex numElem() const
  {
    return m_oldStress.size( 0 );
  }

  /// Number of quadrature points per element storing solid data
  localIndex numQuad() const
  {
    return m_oldStress.size( 1 );
  }

  /**
   * @name Accessors
   */
  ///@{

  /// Non-const/mutable accessor for stress
  arrayView3d< real64, solid::STRESS_USD > const getStress()
  {
    return m_newStress;
  }

  /// Const/non-mutable accessor for stress
  arrayView3d< real64 const, solid::STRESS_USD > const getStress() const
  {
    return m_newStress;
  }

  /// Non-const/Mutable accessor for density.
  arrayView2d< real64 > const getDensity()
  {
    return m_density;
  }

  /// Const/non-mutable accessor for density
  arrayView2d< real64 const > const getDensity() const
  {
    return m_density;
  }

  ///@}

protected:

  virtual void PostProcessInput() override;

  /// The current stress at a quadrature point (i.e. at timestep n, global newton iteration k)
  array3d< real64, solid::STRESS_PERMUTATION > m_newStress;

  /// The previous stress at a quadrature point (i.e. at timestep (n-1))
  array3d< real64, solid::STRESS_PERMUTATION > m_oldStress;

  /// The material density at a quadrature point.
  array2d< real64 > m_density;

  /// The default density for new allocations.
  real64 m_defaultDensity = 0;

  /// band-aid fix...going to have to remove this after we clean up
  /// initialization for constitutive models.
  bool m_postProcessed = false;
};

} // namespace constitutive
} // namespace geosx

#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_ */
