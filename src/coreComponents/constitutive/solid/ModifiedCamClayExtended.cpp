/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file ModifiedCamClayExtended.cpp
 */

#include "ModifiedCamClayExtended.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ModifiedCamClayExtended::ModifiedCamClayExtended( string const & name, 
                                                  Group * const parent ):
  ElasticIsotropic( name, parent )
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ModifiedCamClayExtended, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */

