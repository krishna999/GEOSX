/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "ChomboCoupler.hpp"
#include "hdf5_interface/coupler.hpp"
#include "mesh/FaceManager.hpp"
#include <cstdint>
#include <tuple>

#include <cstdio>

namespace geosx 
{

ChomboCoupler::ChomboCoupler(MPI_Comm const comm, const std::string& outputPath, const std::string& inputPath, MeshLevel& mesh):
  m_comm(comm),
  m_outputPath(outputPath),
  m_inputPath(inputPath),
  m_face_offset(-1),
  m_n_faces_written(-1),
  m_node_offset(-1),
  m_n_nodes_written(-1),
  m_mesh(mesh),
  m_counter(0)
{
  m_mesh.getFaceManager()->RegisterViewWrapper< array1d<real64> >("ChomboPressure");
}

void ChomboCoupler::write(double dt)
{
  ++m_counter;
  FaceManager* faces = m_mesh.getFaceManager();
  NodeManager* nodes = m_mesh.getNodeManager();

  const OrderedVariableOneToManyRelation& face_connectivity = faces->nodeList();
  const localIndex n_faces = face_connectivity.size();

  /* Copy the face connectivity into a contiguous array. */
  std::int64_t* connectivity_array = new std::int64_t[4 * n_faces];
  for (localIndex i = 0; i < n_faces; ++i)
  {
    for (localIndex j = 0; j < 4; ++j)
    {
      connectivity_array[4 * i + j] = face_connectivity[i][j];
    }
  }

  arrayView1d<integer const> const & ruptureState = faces->getReference<integer_array>("ruptureState");
  arrayView1d<integer const> const & ghostRank = faces->getReference<integer_array>(faces->viewKeys.ghostRank);

  bool * faceMask = new bool[n_faces];
  for (localIndex i = 0; i < n_faces; ++i)
  {
    faceMask[i] = (ruptureState[i] > 1) && (ghostRank[i] < 0);
  }

  /* Build the face FieldMap. */
  FieldMap_in face_fields;
  real64 * pressure_ptr = faces->getReference<real64_array>("ChomboPressure").data();
  face_fields["Pressure"] = std::make_tuple(H5T_NATIVE_DOUBLE, 1, pressure_ptr);

  /* Build the node FieldMap. */
  r1_array const& reference_pos = nodes->getReference<r1_array>(nodes->viewKeys.referencePosition);
  localIndex const n_nodes = reference_pos.size();
  R1Tensor const* reference_pos_ptr = reference_pos.data();

  R1Tensor const* displacement_ptr = nodes->getReference<r1_array>(nodes->viewKeys.totalDisplacement).data();            

  FieldMap_in node_fields;
  node_fields["position"] = std::make_tuple(H5T_NATIVE_DOUBLE, 3, reference_pos_ptr);
  node_fields["displacement"] = std::make_tuple(H5T_NATIVE_DOUBLE, 3, displacement_ptr);

  writeBoundaryFile(m_comm, m_outputPath.data(), dt, faceMask,
    m_face_offset, m_n_faces_written, n_faces, connectivity_array, face_fields,
    m_node_offset, m_n_nodes_written, n_nodes,                     node_fields);

  delete[] connectivity_array;
  delete[] faceMask;
}

void ChomboCoupler::read(bool usePressures)
{
  GEOS_LOG_RANK_0("Waiting for file existence: " << m_inputPath);
  waitForFileExistence(m_comm, m_inputPath.data());

  GEOS_LOG_RANK_0("File found: " << m_inputPath);

  if (usePressures)
  {
    GEOS_LOG_RANK_0("Reading pressures...");

    FaceManager* faces = m_mesh.getFaceManager();
    const localIndex n_faces = faces->size();
    const localIndex n_nodes = m_mesh.getNodeManager()->size();

    /* Build the face FieldMap. */
    FieldMap_out face_fields;
    real64 * pressure_ptr = faces->getReference<real64_array>("ChomboPressure").data();
    face_fields["Pressure"] = std::make_tuple(H5T_NATIVE_DOUBLE, 1, pressure_ptr);

    FieldMap_out node_fields;

    readBoundaryFile(m_comm, m_inputPath.data(),
                     m_face_offset, m_n_faces_written, n_faces, face_fields,
                     m_node_offset, m_n_nodes_written, n_nodes, node_fields);
  }

  int rank;
  MPI_Comm_rank(m_comm, &rank);
  if (rank == 0)
  {
    std::remove(m_inputPath.data());
  }
}

} // namespace geosx