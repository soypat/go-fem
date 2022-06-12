package fem

import "gonum.org/v1/gonum/mat"

type Constituter interface {
	Constitutive() mat.Matrix
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
func (m IsotropicMaterial) Constitutive() mat.Matrix {
	G := m.E / (2 + 2*m.Poisson)
	lambda := m.E * m.Poisson / ((1 + m.Poisson) * (1 - 2*m.Poisson))
	return mat.NewDense(6, 6, []float64{
		lambda + 2*G, lambda, lambda, 0, 0, 0,
		lambda, lambda + 2*G, lambda, 0, 0, 0,
		lambda, lambda, lambda + 2*G, 0, 0, 0,
		0, 0, 0, G, 0, 0,
		0, 0, 0, 0, G, 0,
		0, 0, 0, 0, 0, G,
	})
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
func (m TransverselyIsotropicMaterial) Constitutive() mat.Matrix {
	iEx := 1 / m.Ex
	iExy := 1 / m.Exy
	iGxy := 1 / m.Gxy
	cx := -m.PoissonXY * iEx
	cxy := -m.PoissonYZ * iExy
	Gyz := m.Exy / (2*m.PoissonYZ + 2)
	// S is the compliance matrix.
	S := mat.NewDense(6, 6, []float64{
		iEx, cx, cx, 0, 0, 0,
		cx, iExy, cxy, 0, 0, 0,
		cx, cxy, iExy, 0, 0, 0,
		0, 0, 0, 1 / Gyz, 0, 0,
		0, 0, 0, 0, iGxy, 0,
		0, 0, 0, 0, 0, iGxy,
	})
	// Calculate stiffness/constitutive from compliance.
	err := S.Inverse(S)
	if err != nil {
		panic("singular matrix in constitutive calculation probably due to invalid material parameter")
	}
	return S
}
