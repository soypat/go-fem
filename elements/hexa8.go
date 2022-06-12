package elements

import (
	"github.com/soypat/go-fem"

	"gonum.org/v1/gonum/spatial/r3"
)

// Hexa8 is the 3D linear strain hexahedral element of 8 nodes.
type Hexa8 struct{}

var _ fem.Isoparametric3 = Hexa8{}

// LenNodes returns the number of nodes in the element.
func (Hexa8) LenNodes() int { return 8 }

func (Hexa8) Dofs() fem.DofsFlag {
	return fem.DofU
}

// IsoparametricNodes returns the positions of the nodes relative to the origin of the element.
func (Hexa8) IsoparametricNodes() []r3.Vec {
	return []r3.Vec{
		0: {X: -1, Y: -1, Z: -1},
		1: {X: 1, Y: -1, Z: -1},
		2: {X: 1, Y: 1, Z: -1},
		3: {X: -1, Y: 1, Z: -1},
		4: {X: -1, Y: -1, Z: 1},
		5: {X: 1, Y: -1, Z: 1},
		6: {X: 1, Y: 1, Z: 1},
		7: {X: -1, Y: 1, Z: 1},
	}
}

// Basis returns the form functions of the Hexa8 element evaluated at v.
func (Hexa8) Basis(v r3.Vec) []float64 {
	return []float64{
		0: (v.Z*v.Y - v.Y - v.X - v.Z + v.Z*v.X + v.Y*v.X - v.Z*v.Y*v.X + 1) / 8,
		1: (v.X - v.Y - v.Z + v.Z*v.Y - v.Z*v.X - v.Y*v.X + v.Z*v.Y*v.X + 1) / 8,
		2: (v.Y - v.Z + v.X - v.Z*v.Y - v.Z*v.X + v.Y*v.X - v.Z*v.Y*v.X + 1) / 8,
		3: (v.Y - v.Z - v.X - v.Z*v.Y + v.Z*v.X - v.Y*v.X + v.Z*v.Y*v.X + 1) / 8,
		4: (v.Z - v.Y - v.X - v.Z*v.Y - v.Z*v.X + v.Y*v.X + v.Z*v.Y*v.X + 1) / 8,
		5: (v.Z - v.Y + v.X - v.Z*v.Y + v.Z*v.X - v.Y*v.X - v.Z*v.Y*v.X + 1) / 8,
		6: (v.Z + v.Y + v.X + v.Z*v.Y + v.Z*v.X + v.Y*v.X + v.Z*v.Y*v.X + 1) / 8,
		7: (v.Z + v.Y - v.X + v.Z*v.Y - v.Z*v.X - v.Y*v.X - v.Z*v.Y*v.X + 1) / 8,
	}
}

// BasisDiff returns the differentiated form functions of the Hexa8 element
// evaluated at v. The result first contains the form functions differentiated
// with respect to X, then Y and finally Z.
//      [ dN/dx ]
//  N = | dN/dy |   (row major form)
//      [ dN/dz ]
func (Hexa8) BasisDiff(v r3.Vec) []float64 {
	return []float64{
		// Differentiated w.r.t. X
		0: (v.Z + v.Y - (v.Z * v.Y) - 1) / 8,
		1: ((v.Z * v.Y) - v.Y - v.Z + 1) / 8,
		2: (v.Y - v.Z - (v.Z * v.Y) + 1) / 8,
		3: (v.Z - v.Y + (v.Z * v.Y) - 1) / 8,
		4: (v.Y - v.Z + (v.Z * v.Y) - 1) / 8,
		5: (v.Z - v.Y - (v.Z * v.Y) + 1) / 8,
		6: (v.Z + v.Y + (v.Z * v.Y) + 1) / 8,
		7: (-v.Z - v.Y - (v.Z * v.Y) - 1) / 8,
		// Differentiated w.r.t. Y
		8:  (v.Z + v.X - (v.Z * v.X) - 1) / 8,
		9:  (v.Z - v.X + (v.Z * v.X) - 1) / 8,
		10: (v.X - v.Z - (v.Z * v.X) + 1) / 8,
		11: ((v.Z * v.X) - v.X - v.Z + 1) / 8,
		12: (v.X - v.Z + (v.Z * v.X) - 1) / 8,
		13: (-v.Z - v.X - (v.Z * v.X) - 1) / 8,
		14: (v.Z + v.X + (v.Z * v.X) + 1) / 8,
		15: (v.Z - v.X - (v.Z * v.X) + 1) / 8,
		// Differentiated w.r.t. Z
		16: (v.Y + v.X - (v.Y * v.X) - 1) / 8,
		17: (v.Y - v.X + (v.Y * v.X) - 1) / 8,
		18: (-v.Y - v.X - (v.Y * v.X) - 1) / 8,
		19: (v.X - v.Y + (v.Y * v.X) - 1) / 8,
		20: ((v.Y * v.X) - v.X - v.Y + 1) / 8,
		21: (v.X - v.Y - (v.Y * v.X) + 1) / 8,
		22: (v.Y + v.X + (v.Y * v.X) + 1) / 8,
		23: (v.Y - v.X - (v.Y * v.X) + 1) / 8,
	}
}

func (Hexa8) Quadrature() (positions []r3.Vec, weights []float64) {
	positions, weights, _ = uniformGaussQuad(2, 2, 2)
	return positions, weights
}
