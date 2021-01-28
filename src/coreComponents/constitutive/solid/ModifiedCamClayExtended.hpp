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
 *  @file ModifiedCamClayExtended.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_MODIFIEDCAMCLAYEXTENDED_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_MODIFIEDCAMCLAYEXTENDED_HPP_

#include "ElasticIsotropic.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class ModifiedCamClayExtendedUpdates
 *
 * Class to provide material updates that may be called from a kernel function.
 */

class ModifiedCamClayExtendedUpdates : public ElasticIsotropicUpdates
{
public:

  ModifiedCamClayExtendedUpdates( arrayView1d< real64 const > const & bulkModulus,
                                  arrayView1d< real64 const > const & shearModulus,
                                  arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                  arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    ElasticIsotropicUpdates( bulkModulus, 
                             shearModulus, 
                             newStress, 
                             oldStress )
  {}

  /// Default copy constructor

  ModifiedCamClayExtendedUpdates( ModifiedCamClayExtendedUpdates const & ) = default;

  /// Default move constructor

  ModifiedCamClayExtendedUpdates( ModifiedCamClayExtendedUpdates && ) = default;

  /// Deleted default constructor

  ModifiedCamClayExtendedUpdates() = delete;

  /// Deleted copy assignment operator

  ModifiedCamClayExtendedUpdates & operator = ( ModifiedCamClayExtendedUpdates const & ) = delete;

  /// Deleted move assignment operator

  ModifiedCamClayExtendedUpdates & operator = ( ModifiedCamClayExtendedUpdates && ) = delete;

};

/**
 * @class ModifiedCamClayExtended
 *
 * Extended modified Cam-Clay model for hard rock.
 */

class ModifiedCamClayExtended : public ElasticIsotropic
{
public:

  ModifiedCamClayExtended( string const & name, 
                           Group * const parent );

  /**
   * Default Destructor
   */

  virtual ~ModifiedCamClayExtended() override = default;

  /// string name to use for this class in the catalog

  static constexpr auto m_catalogNameString = "ModifiedCamClayExtended";

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_MODIFIEDCAMCLAYEXTENDED_HPP_ */

