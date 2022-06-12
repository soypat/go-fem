package elements

import (
	"github.com/soypat/go-fem"

	"gonum.org/v1/gonum/spatial/r3"
)

// Tetra10 is the 3D quadratic strain tetrahedral element of 4 corner nodes
// and 6 edge nodes.
type Tetra10 struct{}

var _ fem.Isoparametric3 = Tetra10{}

func (Tetra10) Dofs() fem.DofsFlag {
	return fem.DofU
}

func (Tetra10) LenNodes() int { return 10 }

// IsoparametricNodes returns the positions of the nodes relative to the origin of the element.
func (Tetra10) IsoparametricNodes() []r3.Vec {
	return []r3.Vec{
		0: {X: 0, Y: 0, Z: 0},
		1: {X: 1, Y: 0, Z: 0},
		2: {X: 0, Y: 1, Z: 0},
		3: {X: 0, Y: 0, Z: 1},
		4: {X: 1 / 2, Y: 0, Z: 0},
		5: {X: 0, Y: 1 / 2, Z: 0},
		6: {X: 0, Y: 0, Z: 1 / 2},
		7: {X: 1 / 2, Y: 1 / 2, Z: 0},
		8: {X: 0, Y: 1 / 2, Z: 1 / 2},
		9: {X: 1 / 2, Y: 0, Z: 1 / 2},
	}
}

// Basis returns the form functions of the Tetra10 element evaluated at v.
func (Tetra10) Basis(v r3.Vec) []float64 {
	return []float64{
		0: 2*v.X*v.X + 4*v.X*v.Y + 4*v.X*v.Z - 3*v.X + 2*v.Y*v.Y + 4*v.Y*v.Z - 3*v.Y + 2*v.Z*v.Z - 3*v.Z + 1,
		1: 2*v.X*v.X - v.X,
		2: 2*v.Y*v.Y - v.Y,
		3: 2*v.Z*v.Z - v.Z,
		4: 4*v.X - 4*v.X*v.Y - 4*v.X*v.Z - 4*v.X*v.X,
		5: 4*v.Y - 4*v.X*v.Y - 4*v.Y*v.Z - 4*v.Y*v.Y,
		6: 4*v.Z - 4*v.X*v.Z - 4*v.Y*v.Z - 4*v.Z*v.Z,
		7: 4 * v.X * v.Y,
		8: 4 * v.Y * v.Z,
		9: 4 * v.X * v.Z,
	}
}

func (Tetra10) BasisDiff(v r3.Vec) []float64 {
	return []float64{
		// w.r.t. X
		0: 4*v.X + 4*v.Y + 4*v.Z - 3,
		1: 4*v.X - 1,
		2: 0,
		3: 0,
		4: 4 - 4*v.Y - 4*v.Z - 8*v.X,
		5: -4 * v.Y,
		6: -4 * v.Z,
		7: 4 * v.Y,
		8: 0,
		9: 4 * v.Z,
		// w.r.t. Y
		10: 4*v.X + 4*v.Y + 4*v.Z - 3,
		11: 0,
		12: 4*v.Y - 1,
		13: 0,
		14: -4 * v.X,
		15: 4 - 8*v.Y - 4*v.Z - 4*v.X,
		16: -4 * v.Z,
		17: 4 * v.X,
		18: 4 * v.Z,
		19: 0,
		// w.r.t. Z
		20: 4*v.X + 4*v.Y + 4*v.Z - 3,
		21: 0,
		22: 0,
		23: 4*v.Z - 1,
		24: -4 * v.X,
		25: -4 * v.Y,
		26: 4 - 4*v.Y - 8*v.Z - 4*v.X,
		27: 0,
		28: 4 * v.Y,
		29: 4 * v.X,
	}
}

func (Tetra10) Quadrature() (positions []r3.Vec, weights []float64) {
	const (
		sqrt5 = 2.2360679774997896964091736687312762354
		a     = (5 + 3*sqrt5) / 20
		b     = (5 - sqrt5) / 20
	)
	positions = []r3.Vec{
		{X: a, Y: b, Z: b},
		{X: b, Y: b, Z: b},
		{X: b, Y: b, Z: a},
		{X: b, Y: a, Z: b},
	}
	weights = []float64{0.25, 0.25, 0.25, 0.25}
	return positions, weights
}
