package fem

import (
	"fmt"

	"gonum.org/v1/gonum/floats/scalar"
	"gonum.org/v1/gonum/mat"
)

// Constituter represents the homogenous properties of a medium
// that can then be used to model solids or other continuous field problems.
// For solids it returns the unmodified constitutive tensor (Generalized Hookes law).
// The shear modulii should not be halved.
type Constituter interface {
	Constitutive() (mat.Matrix, error)
}

// IsotropicMaterial represents a material with no preferential direction
// respecting mechanical properties.
// i.e: Room temperature steel. 200GPa, Poisson=0.3
type IsotropicMaterial struct {
	// Elasiticity Modulus. Also known as Young's modulus.
	E float64
	// Poisson modulus. Usually uses nu as symbol in literature.
	Poisson float64
}

// return value will be concrete in future.
func (m IsotropicMaterial) Constitutive() (mat.Matrix, error) {
	G := m.ShearModulus()
	lambda := m.lameParam2()
	return mat.NewDense(6, 6, []float64{
		lambda + 2*G, lambda, lambda, 0, 0, 0,
		lambda, lambda + 2*G, lambda, 0, 0, 0,
		lambda, lambda, lambda + 2*G, 0, 0, 0,
		0, 0, 0, G, 0, 0,
		0, 0, 0, 0, G, 0,
		0, 0, 0, 0, 0, G,
	}), nil
}

// ShearModulus returns the shear modulus of the material, often denoted as G.
func (m IsotropicMaterial) ShearModulus() float64 {
	return m.E / (2 + 2*m.Poisson)
}

// lameParam2 returns the second Lamé parameter, also known as lambda.
// https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
func (m IsotropicMaterial) lameParam2() float64 {
	return m.E * m.Poisson / ((1 + m.Poisson) * (1 - 2*m.Poisson))
}

// TransverselyIsotropicMaterial is a material that is transversally isotropic
// with X as the principal or longitudinal direction.
// i.e. AS4 carbon fiber: Ex=235GPa, Exy=14GPa, Gxy=28GPa, vxy=0.2, vyz=0.25
type TransverselyIsotropicMaterial struct {
	// Longitudinal Young's elasticity modulus (X direction). Found as E_1 in literature.
	Ex float64
	// Transversal Young's elasticity modulus (XY plane). Found as E_2 in literature.
	Exy float64
	// Transversal shear modulus in XY plane. Is equal to Gxz for this material.
	Gxy float64
	// PoissonXY represents the quotient between deformations in directions X and Y
	// when stress is applied in direction X. Is positive when material contracts in transversal direction.
	PoissonXY float64
	// PoissonXY represents the quotient between deformations in directions Y and Z
	// when stress is applied in direction Y. Is positive when material contracts in transversal direction.
	PoissonYZ float64
}

// return value will be concrete in future.
func (m TransverselyIsotropicMaterial) Constitutive() (mat.Matrix, error) {
	// # Matlab program to solve for inverse of compliance matrix:
	// syms Ex Exy Gxy vxy vyz
	// cx = -vxy/Ex;
	// cxy = -vyz/Exy;
	// Gyz = Exy / (2*vyz + 2);
	// comp = [1/Ex, cx, cx, 0, 0, 0
	// 		cx, 1/Exy, cxy, 0, 0, 0
	// 		cx, cxy, 1/Exy, 0, 0, 0
	// 		0, 0, 0, 1 / Gyz, 0, 0
	// 		0, 0, 0, 0, 1/ Gxy, 0
	// 		0, 0, 0, 0, 0, 1 / Gxy,];
	// inv(comp)
	const tol = 1e-10
	Ex := m.Ex
	Exy := m.Exy
	Gxy := m.Gxy
	vxy := m.PoissonXY
	vyz := m.PoissonYZ
	vxy2 := vxy * vxy
	vyz2 := vyz * vyz
	div1 := (2*Exy*vxy2 - Ex + Ex*vyz)
	div2 := (2*Exy*vxy2*vyz + 2*Exy*vxy2 + Ex*vyz2 - Ex)
	if scalar.EqualWithinAbs(div1, 0, tol) || scalar.EqualWithinAbs(div2, 0, tol) {
		return nil, fmt.Errorf("material parameters yielded singular matrix")
	}
	data := []float64{
		(Ex*Ex*vyz - Ex*Ex) / div1, -(Ex * Exy * vxy) / div1, -(Ex * Exy * vxy) / div1, 0, 0, 0,
		-(Ex * Exy * vxy) / div1, -(Exy * (-Exy*vxy2 + Ex)) / div2, -(Exy * (Exy*vxy2 + Ex*vyz)) / div2, 0, 0, 0,
		-(Ex * Exy * vxy) / div1, -(Exy * (Exy*vxy2 + Ex*vyz)) / div2, -(Exy * (-Exy*vxy2 + Ex)) / div2, 0, 0, 0,
		0, 0, 0, Exy / (2 * (vyz + 1)), 0, 0,
		0, 0, 0, 0, Gxy, 0,
		0, 0, 0, 0, 0, Gxy,
	}
	return mat.NewDense(6, 6, data), nil
}

type AxisymmetricIsotropicMaterial struct {
	IsotropicMaterial
}

func (ax AxisymmetricIsotropicMaterial) Constitutive() (mat.Matrix, error) {
	nu := ax.Poisson
	f := nu / (1 - nu)
	g := (1 - 2*nu) / (2 * (1 - nu))
	factor := ax.E * (1 - nu) / ((1 + nu) * (1 - 2*nu))
	f *= factor
	g *= factor
	data := []float64{
		factor, f, f, 0, 0, 0,
		f, factor, f, 0, 0, 0,
		f, f, factor, 0, 0, 0,
		0, 0, 0, g, 0, 0,
		0, 0, 0, 0, g, 0,
		0, 0, 0, 0, 0, g}
	d := mat.NewDense(6, 6, data)
	return d, nil
}
