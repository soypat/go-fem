package solids

import (
	"math"

	"gonum.org/v1/gonum/floats/scalar"
	"gonum.org/v1/gonum/mat"
)

// TransverselyIsotropic is a material that is transversally isotropic
// with X as the principal or longitudinal direction.
// i.e. AS4 carbon fiber: Ex=235GPa, Exy=14GPa, Gxy=28GPa, vxy=0.2, vyz=0.25
type TransverselyIsotropic struct {
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
func (m TransverselyIsotropic) Constitutive() (mat.Matrix, error) {
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
		return nil, mat.Condition(math.Min(div1, div2))
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
