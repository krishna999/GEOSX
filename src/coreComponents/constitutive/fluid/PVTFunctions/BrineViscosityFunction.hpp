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
 * @file BrineViscosityFunction.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEVISCOSITYFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEVISCOSITYFUNCTION_HPP_

#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class BrineViscosityFunction : public PVTFunction
{
public:


  BrineViscosityFunction( string_array const & inputPara,
                          string_array const & componentNames,
                          real64_array const & componentMolarWeight );
  ~BrineViscosityFunction() override
  {}

  static constexpr auto m_catalogName = "BrineViscosity";
  static string catalogName()                    { return m_catalogName; }
  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTFuncType functionType() const override
  {
    return PVTFuncType::VISCOSITY;

  }

  virtual void evaluation( EvalVarArgs const & pressure,
                           EvalVarArgs const & temperature,
                           arraySlice1d< EvalVarArgs const > const & phaseComposition,
                           EvalVarArgs & value, bool useMass = 0 ) const override;

private:

  void makeCoef( string_array const & inputPara );

  real64 m_coef0;
  real64 m_coef1;

};

}

}
#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEVISCOSITYFUNCTION_HPP_
