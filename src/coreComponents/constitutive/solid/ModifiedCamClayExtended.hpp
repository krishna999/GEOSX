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
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class ModifiedCamClayExtendedUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class ModifiedCamClayExtendedUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] friction The ArrayView holding the friction data for each element.
   * @param[in] dilation The ArrayView holding the dilation data for each element.
   * @param[in] hardening The ArrayView holding the hardening data for each element.
   * @param[in] newCohesion The ArrayView holding the new cohesion data for each element.
   * @param[in] oldCohesion The ArrayView holding the old cohesion data for each element.
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  ModifiedCamClayExtendedUpdates( arrayView1d< real64 const > const & friction,
                        arrayView1d< real64 const > const & dilation,
                        arrayView1d< real64 const > const & hardening,
                        arrayView2d< real64 > const & newCohesion,
                        arrayView2d< real64 > const & oldCohesion,
                        arrayView1d< real64 const > const & bulkModulus,
                        arrayView1d< real64 const > const & shearModulus,
                        arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    ElasticIsotropicUpdates( bulkModulus, shearModulus, newStress, oldStress ),
    m_friction( friction ),
    m_dilation( dilation ),
    m_hardening( hardening ),
    m_newCohesion( newCohesion ),
    m_oldCohesion( oldCohesion )
  {}

  /// Default copy constructor
  ModifiedCamClayExtendedUpdates( ModifiedCamClayExtendedUpdates const & ) = default;

  /// Default move constructor
  ModifiedCamClayExtendedUpdates( ModifiedCamClayExtendedUpdates && ) = default;

  /// Deleted default constructor
  ModifiedCamClayExtendedUpdates() = delete;

  /// Deleted copy assignment operator
  ModifiedCamClayExtendedUpdates & operator=( ModifiedCamClayExtendedUpdates const & ) = delete;

  /// Deleted move assignment operator
  ModifiedCamClayExtendedUpdates & operator=( ModifiedCamClayExtendedUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic;

  // Bring in base implementations to prevent hiding warnings
  using ElasticIsotropicUpdates::smallStrainUpdate;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  real64 ( &stiffness )[6][6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const final;

private:
  /// A reference to the ArrayView holding the friction angle for each element.
  arrayView1d< real64 const > const m_friction;

  /// A reference to the ArrayView holding the dilation angle for each element.
  arrayView1d< real64 const > const m_dilation;

  /// A reference to the ArrayView holding the hardening rate for each element.
  arrayView1d< real64 const > const m_hardening;

  /// A reference to the ArrayView holding the new cohesion for each integration point
  arrayView2d< real64 > const m_newCohesion;

  /// A reference to the ArrayView holding the old cohesion for each integration point
  arrayView2d< real64 > const m_oldCohesion;


  // TODO to change to MMC yield function with cohesion

  GEOSX_HOST_DEVICE
  void yieldAndDerivatives( real64 const invP,
                            real64 const invQ,
                            real64 const friction,
                            real64 const cohesion,
                            real64 (& FdF)[4] ) const
  {
    real64 F     =  invQ + friction * invP - cohesion;
    real64 dF_dP = friction;
    real64 dF_dQ = 1.0;
    real64 dF_dCohesion = -1.0; 

    FdF[0] = F;
    FdF[1] = dF_dP;
    FdF[2] = dF_dQ;
    FdF[3] = dF_dCohesion;
  }

  // TODO to change to MMC hardening function

  GEOSX_HOST_DEVICE
  void hardeningAndDerivatives( real64 const cohesion,
                                real64 const plasticMultiplier,
                                real64 const hardeningRate,
                                real64 (& HdH)[2] ) const
  {
    real64 H = cohesion + plasticMultiplier * hardeningRate;
    real64 dH_dPlasticMultiplier = hardeningRate; 

    if( H < 0 )
    {
      H = 0;
      dH_dPlasticMultiplier = 0;
    }

    HdH[0] = H;
    HdH[1] = dH_dPlasticMultiplier;
  }

  // TODO to change to MMC plastic potential function

  GEOSX_HOST_DEVICE
  void potentialDerivatives( real64 const dilation,
                             real64 (& dG)[2] ) const
  {
    // G = Q + dilation * P
    real64 dG_dP = dilation;
    real64 dG_dQ = 1.0; 

    dG[0] = dG_dP;
    dG[1] = dG_dQ;
  }

};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ModifiedCamClayExtendedUpdates::smallStrainUpdate( localIndex const k,
                                              localIndex const q,
                                              real64 const ( &strainIncrement )[6],
                                              real64 ( & stress )[6],
                                              real64 ( & stiffness )[6][6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)

  ElasticIsotropicUpdates::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

  // decompose into mean (P) and von mises (Q) stress invariants

  real64 trialP;
  real64 trialQ;
  real64 deviator[6];

  twoInvariant::stressDecomposition( stress,
                                     trialP,
                                     trialQ,
                                     deviator );

  // check yield function F <= 0, using old hardening variable state
  real64 FdF[4];
  yieldAndDerivatives( trialP,
                       trialQ,
                       m_friction[k],
                       m_oldCohesion[k][q],
                       FdF );

  real64 yield = FdF[0];//trialQ + m_friction[k] * trialP - m_oldCohesion[k][q];

  if( yield < 1e-9 ) // elasticity
  {
    return;
  }

  // else, plasticity (trial stress point lies outside yield surface)

  // the return mapping can in general be written as a newton iteration.
  // here we have a linear problem, so the algorithm will converge in one
  // iteration, but this is a template for more general models with either
  // nonlinear hardening or yield surfaces.

  real64 solution[3], residual[3], delta[3];
  real64 jacobian[3][3] = {{}}, jacobianInv[3][3] = {{}};

  solution[0] = trialP; // initial guess for newP
  solution[1] = trialQ; // initial guess for newQ
  solution[2] = 1e-5;   // initial guess for plastic multiplier

  real64 norm, normZero = 1e30;

  // begin newton loop

  for( localIndex iter=0; iter<20; ++iter )
  {
    // apply a linear cohesion decay model,
    // then check for complete cohesion loss
    
    real64 HdH[2];
    hardeningAndDerivatives( m_oldCohesion[k][q],
                             solution[2],
                             m_hardening[k],
                             HdH );


    m_newCohesion[k][q]  = HdH[0];//m_oldCohesion[k][q] + solution[2] * m_hardening[k];
    //real64 cohesionDeriv = HdH[1];//m_hardening[k];
/**
    if( m_newCohesion[k][q] < 0 )
    {
      m_newCohesion[k][q] = 0;
      cohesionDeriv = 0;
    }
*/
    // assemble residual system
    // resid1 = P - trialP + dlambda*bulkMod*dG/dP = 0
    // resid2 = Q - trialQ + dlambda*3*shearMod*dG/dQ = 0
    // resid3 = F = 0

    yieldAndDerivatives( solution[0],
                         solution[1],
                         m_friction[k],
                         m_newCohesion[k][q],
                         FdF );

    real64 dG[2];
    potentialDerivatives( m_dilation[k],
                          dG );

    residual[0] = solution[0] - trialP + m_bulkModulus[k] * solution[2] * dG[0];//solution[0] - trialP + solution[2] * m_bulkModulus[k] * m_dilation[k];
    residual[1] = solution[1] - trialQ + 3 * m_shearModulus[k] * solution[2] * dG[1];//solution[1] - trialQ + 3 * solution[2] * m_shearModulus[k];
    residual[2] = FdF[0];//solution[1] + m_friction[k] * solution[0] - m_newCohesion[k][q];

    // check for convergence

    norm = LvArray::tensorOps::l2Norm< 3 >( residual );

    if( iter==0 )
    {
      normZero = norm;
    }

    if( norm < 1e-8*(normZero+1))
    {
      break;
    }

    // solve Newton system

    jacobian[0][0] = 1.0;
    jacobian[0][2] = m_bulkModulus[k] * dG[0];
    jacobian[1][1] = 1.0;
    jacobian[1][2] = 3.0 * m_shearModulus[k] * dG[1];
    jacobian[2][0] = FdF[1];//m_friction[k];
    jacobian[2][1] = FdF[2];//1;
    jacobian[2][2] = FdF[3] * HdH[1];//cohesionDeriv;

    LvArray::tensorOps::invert< 3 >( jacobianInv, jacobian );
    LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( delta, jacobianInv, residual );

    for( localIndex i=0; i<3; ++i )
    {
      solution[i] -= delta[i];
    }
  }

  // re-construct stress = P*eye + sqrt(2/3)*Q*nhat

  twoInvariant::stressRecomposition( solution[0],
                                     solution[1],
                                     deviator,
                                     stress );

  // construct consistent tangent stiffness
  // note: if trialQ = 0, we will get a divide by zero error below,
  // but this is an unphysical (zero-strength) state anyway

  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );

  real64 c1 = 2 * m_shearModulus[k] * solution[1] / trialQ;
  real64 c2 = jacobianInv[0][0] * m_bulkModulus[k] - c1 / 3;
  real64 c3 = sqrt( 2./3 ) * 3 * m_shearModulus[k] * jacobianInv[0][1];
  real64 c4 = sqrt( 2./3 ) * m_bulkModulus[k] * jacobianInv[1][0];
  real64 c5 = 2 * jacobianInv[1][1] * m_shearModulus[k] - c1;

  real64 identity[6];

  for( localIndex i=0; i<3; ++i )
  {
    stiffness[i][i] = c1;
    stiffness[i+3][i+3] = 0.5 * c1;
    identity[i] = 1.0;
    identity[i+3] = 0.0;
  }

  for( localIndex i=0; i<6; ++i )
  {
    for( localIndex j=0; j<6; ++j )
    {
      stiffness[i][j] +=   c2 * identity[i] * identity[j]
                         + c3 * identity[i] * deviator[j]
                         + c4 * deviator[i] * identity[j]
                         + c5 * deviator[i] * deviator[j];
    }
  }

  if(k==0)
  {
    

    // comparison with zhang 1994, appendixes A, B
    real64 PP = FdF[2];
    real64 QQ = FdF[1];
    real64 A11 = PP;
    real64 A12 = QQ;
    real64 A21 = 0.;
    real64 A22 = 0.;

    real64 B11 = 0.;
    real64 B12 = 0.;
    real64 B21 = QQ/3.;
    real64 B22 = -PP;

    real64 KK = m_bulkModulus[k];
    real64 GG = m_shearModulus[k];

    real64 Delta = (A11 + 3. * KK * B11) * (A22 + 3. * GG * B22) - (A12 + 3. * GG * B12) * (A21 + 3. * KK * B21);

    real64 C11 = ((A22 + 3. * GG * B22)*B11 - (A12 + 3. * GG * B12)*B21) / Delta;
    real64 C21 = ((A11 + 3. * KK * B11)*B21 - (A21 + 3. * KK * B21)*B11) / Delta;
    real64 C12 = ((A22 + 3. * GG * B22)*B12 - (A12 + 3. * GG * B12)*B22) / Delta;
    real64 C22 = ((A11 + 3. * KK * B11)*B22 - (A21 + 3. * KK * B21)*B12) / Delta;

    real64 d1 = c1;
    real64 d2 = KK - 3.0 * KK * KK * C11 - d1/3.0;  
    real64 d3 = -sqrt(3.0/2.0) * ( -2.*GG*KK*C12 );  // factor -sqrt(3/2) is because of the difference in normalized deviatoric
    real64 d4 = -sqrt(3.0/2.0) * ( -6.*GG*KK*C21 );  // factor -sqrt(3/2) is because of the difference in normalized deviatoric
    real64 d5 = 3.0/2.0 * ( 4.0/3.0 * GG - 2.0/3.0 * d1 - 4.0 * GG * GG * C22 ); // factor 3/2 
                                                                                   //is because of the difference 
                                                                                   //in normalized deviatoric

    // The validation was perfect
    std::cout << "d2/c2 = " << d2/c2 << std::endl;
    std::cout << "d3/c3 = " << d3/c3 << std::endl;
    std::cout << "d4/c4 = " << d4/c4 << std::endl;
    std::cout << "d5/c5 = " << d5/c5 << std::endl;


  }
  // save new stress and return
  saveStress( k, q, stress );
  return;
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ModifiedCamClayExtendedUpdates::smallStrainUpdate( localIndex const k,
                                              localIndex const q,
                                              real64 const ( &strainIncrement )[6],
                                              real64 ( & stress )[6],
                                              DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, strainIncrement, stress, stiffness.m_c );
}



/**
 * @class ModifiedCamClayExtended
 *
 * Modified Cam-Clay extended material model.
 */
class ModifiedCamClayExtended : public ElasticIsotropic
{
public:

  /// @typedef Alias for ModifiedCamClayExtendedUpdates
  using KernelWrapper = ModifiedCamClayExtendedUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  ModifiedCamClayExtended( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~ModifiedCamClayExtended() override;


  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ModifiedCamClayExtended";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default friction angle
    static constexpr auto defaultFrictionAngleString = "defaultFrictionAngle";

    /// string/key for default dilation angle
    static constexpr auto defaultDilationAngleString = "defaultDilationAngle";

    /// string/key for default hardening rate
    static constexpr auto defaultHardeningString = "defaultHardeningRate";

    /// string/key for default cohesion
    static constexpr auto defaultCohesionString = "defaultCohesion";

    /// string/key for friction angle
    static constexpr auto frictionString  = "friction";

    /// string/key for dilation angle
    static constexpr auto dilationString  = "dilation";

    /// string/key for cohesion
    static constexpr auto hardeningString  = "hardening";

    /// string/key for cohesion
    static constexpr auto newCohesionString  = "cohesion";

    /// string/key for cohesion
    static constexpr auto oldCohesionString  = "oldCohesion";
  };

  /**
   * @brief Create a instantiation of the ModifiedCamClayExtendedUpdate class that refers to the data in this.
   * @return An instantiation of ModifiedCamClayExtendedUpdate.
   */
  ModifiedCamClayExtendedUpdates createKernelUpdates() const
  {
    return ModifiedCamClayExtendedUpdates( m_friction,
                                 m_dilation,
                                 m_hardening,
                                 m_newCohesion,
                                 m_oldCohesion,
                                 m_bulkModulus,
                                 m_shearModulus,
                                 m_newStress,
                                 m_oldStress );
  }

protected:
  virtual void postProcessInput() override;

  /// Material parameter: The default value of yield surface slope
  real64 m_defaultFrictionAngle;

  /// Material parameter: The default value of plastic potential slope
  real64 m_defaultDilationAngle;

  /// Material parameter: The default value of the initial cohesion
  real64 m_defaultCohesion;

  /// Material parameter: The default value of the hardening rate
  real64 m_defaultHardening;

  /// Material parameter: The yield surface slope for each element
  array1d< real64 > m_friction;

  /// Material parameter: The plastic potential slope for each element
  array1d< real64 > m_dilation;

  /// Material parameter: The hardening rate each element
  array1d< real64 > m_hardening;

  /// State variable: The current cohesion parameter for each quadrature point
  array2d< real64 > m_newCohesion;

  /// State variable: The previous cohesion parameter for each quadrature point
  array2d< real64 > m_oldCohesion;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_MODIFIEDCAMCLAYEXTENDED_HPP_ */

