package solids

import (
	"math"

	"github.com/soypat/go-fem"
	"gonum.org/v1/gonum/mat"
)

// Isotropic represents a material with no preferential direction
// respecting mechanical properties.
// i.e: Room temperature steel. E=200GPa, Poisson=0.3
type Isotropic struct {
	// Elasiticity Modulus. Also known as Young's modulus.
	E float64
	// Poisson modulus. Usually uses nu as symbol in literature.
	Poisson float64
}

// return value will be concrete in future.
func (m Isotropic) Constitutive() (mat.Matrix, error) {
	G := m.ShearModulus()
	lambda := m.lameParam2()
	if math.IsNaN(lambda) || math.IsInf(lambda, 0) {
		return nil, mat.Condition((1 + m.Poisson) * (1 - 2*m.Poisson))
	}
	return mat.NewDense(6, 6, []float64{
		lambda + 2*G, lambda, lambda, 0, 0, 0,
		lambda, lambda + 2*G, lambda, 0, 0, 0,
		lambda, lambda, lambda + 2*G, 0, 0, 0,
		0, 0, 0, G, 0, 0,
		0, 0, 0, 0, G, 0,
		0, 0, 0, 0, 0, G,
	}), nil
}

func (m Isotropic) Solid3D() fem.IsoConstituter {
	var isoc isoconstituter
	isoc.C, isoc.err = m.Constitutive()
	isoc.strain = SetStrainDisplacementMatrixXYZ
	return isoc
}

// ShearModulus returns the shear modulus of the material, often denoted as G.
func (m Isotropic) ShearModulus() float64 {
	return m.E / (2 + 2*m.Poisson)
}

// lameParam2 returns the second Lamé parameter, also known as lambda.
// https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
func (m Isotropic) lameParam2() float64 {
	return m.E * m.Poisson / ((1 + m.Poisson) * (1 - 2*m.Poisson))
}

func (m Isotropic) PlaneStess() fem.IsoConstituter {
	E := m.E
	nu := m.Poisson
	factor := E / (1 - nu*nu)
	return isoconstituter{
		C: mat.NewDense(3, 3, []float64{
			factor, factor * nu, 0,
			factor * nu, factor, 0,
			0, 0, factor * (1 - nu) / 2,
		}),
		strain: SetStrainDisplacementMatrixPlane,
	}
}

func (m Isotropic) PlaneStrain() fem.IsoConstituter {
	E := m.E
	nu := m.Poisson
	factor := E / (1 + nu) / (1 - 2*nu)
	nuf := nu * factor
	return isoconstituter{
		C: mat.NewDense(3, 3, []float64{
			factor * (1 - nu), nuf, nuf,
			nuf, factor * (1 - nu), nuf,
			nuf, nuf, factor * (1 - 2*nu),
		}),
		strain: SetStrainDisplacementMatrixPlane,
	}
}

func (m Isotropic) Axisymmetric() fem.IsoConstituter {
	nu := m.Poisson
	f := nu / (1 - nu)
	g := (1 - 2*nu) / (2 * (1 - nu))
	factor := m.E * (1 - nu) / ((1 + nu) * (1 - 2*nu))
	f *= factor
	g *= factor
	data := []float64{
		factor, f, f, 0,
		f, factor, f, 0,
		f, f, factor, 0,
		0, 0, 0, g,
	}
	d := mat.NewDense(4, 4, data)
	return isoconstituter{
		C:      d,
		strain: SetStrainDisplacementMatrixAxisymmetric,
	}
}
