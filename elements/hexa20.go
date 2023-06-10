package elements

import (
	"github.com/soypat/go-fem"

	"gonum.org/v1/gonum/spatial/r3"
)

// Hexa20 is the serendipity 3D quadratic strain hexahedral element of 20 nodes.
type Hexa20 struct{}

var _ fem.Isoparametric3 = Hexa20{}

// LenNodes returns the number of nodes in the element.
func (Hexa20) LenNodes() int { return 20 }

// Dofs returns the set dofs for the element's nodes.
func (Hexa20) Dofs() fem.DofsFlag {
	return fem.DofPos
}

// IsoparametricNodes returns the positions of the nodes relative to the origin of the element.
func (Hexa20) IsoparametricNodes() []r3.Vec {
	return []r3.Vec{
		0:  {X: -1, Y: -1, Z: -1},
		1:  {X: -1, Y: 1, Z: -1},
		2:  {X: 1, Y: 1, Z: -1},
		3:  {X: 1, Y: -1, Z: -1},
		4:  {X: -1, Y: -1, Z: 1},
		5:  {X: -1, Y: 1, Z: 1},
		6:  {X: 1, Y: 1, Z: 1},
		7:  {X: 1, Y: -1, Z: 1},
		8:  {X: -1, Y: 0, Z: -1},
		9:  {X: 0, Y: 1, Z: -1},
		10: {X: 1, Y: 0, Z: -1},
		11: {X: 0, Y: -1, Z: -1},
		12: {X: -1, Y: 0, Z: 1},
		13: {X: 0, Y: 1, Z: 1},
		14: {X: 1, Y: 0, Z: 1},
		15: {X: 0, Y: -1, Z: 1},
		16: {X: -1, Y: -1, Z: 0},
		17: {X: -1, Y: 1, Z: 0},
		18: {X: 1, Y: 1, Z: 0},
		19: {X: 1, Y: -1, Z: 0},
	}
}

// Basis returns the form functions of the Hexa20 element evaluated at v.
func (Hexa20) Basis(v r3.Vec) []float64 {
	x2 := v.X * v.X
	y2 := v.Y * v.Y
	z2 := v.Z * v.Z
	return []float64{
		((v.Y - 1) * (v.X - 1) * (v.Z - 1) * (v.Y + v.X + v.Z + 2)) / 8,
		-((v.Y + 1) * (v.X - 1) * (v.Z - 1) * (v.X - v.Y + v.Z + 2)) / 8,
		-((v.Y + 1) * (v.X + 1) * (v.Z - 1) * (v.Y + v.X - v.Z - 2)) / 8,
		-((v.Y - 1) * (v.X + 1) * (v.Z - 1) * (v.Y - v.X + v.Z + 2)) / 8,
		-((v.Y - 1) * (v.X - 1) * (v.Z + 1) * (v.Y + v.X - v.Z + 2)) / 8,
		-((v.Y + 1) * (v.X - 1) * (v.Z + 1) * (v.Y - v.X + v.Z - 2)) / 8,
		((v.Y + 1) * (v.X + 1) * (v.Z + 1) * (v.Y + v.X + v.Z - 2)) / 8,
		((v.Y - 1) * (v.X + 1) * (v.Z + 1) * (v.Y - v.X - v.Z + 2)) / 8,
		-((y2 - 1) * (v.X - 1) * (v.Z - 1)) / 4,
		((x2 - 1) * (v.Y + 1) * (v.Z - 1)) / 4,
		((y2 - 1) * (v.X + 1) * (v.Z - 1)) / 4,
		-((x2 - 1) * (v.Y - 1) * (v.Z - 1)) / 4,
		((y2 - 1) * (v.X - 1) * (v.Z + 1)) / 4,
		-((x2 - 1) * (v.Y + 1) * (v.Z + 1)) / 4,
		-((y2 - 1) * (v.X + 1) * (v.Z + 1)) / 4,
		((x2 - 1) * (v.Y - 1) * (v.Z + 1)) / 4,
		-((z2 - 1) * (v.Y - 1) * (v.X - 1)) / 4,
		((z2 - 1) * (v.Y + 1) * (v.X - 1)) / 4,
		-((z2 - 1) * (v.Y + 1) * (v.X + 1)) / 4,
		((z2 - 1) * (v.Y - 1) * (v.X + 1)) / 4,
	}
}

// BasisDiff returns the differentiated form functions of the Hexa20 element
// evaluated at v. The result first contains the form functions differentiated
// with respect to X, then Y and finally Z.
//
//	    [ dN/dx ]
//	N = | dN/dy |   (row major form)
//	    [ dN/dz ]
func (Hexa20) BasisDiff(v r3.Vec) []float64 {
	x2 := v.X * v.X
	y2 := v.Y * v.Y
	z2 := v.Z * v.Z
	return []float64{
		// Differentiated w.r.t X.
		0:  ((v.Y - 1) * (v.Z - 1) * (v.Y + 2*v.X + v.Z + 1)) / 8,
		1:  -((v.Y + 1) * (v.Z - 1) * (2*v.X - v.Y + v.Z + 1)) / 8,
		2:  -((v.Y + 1) * (v.Z - 1) * (v.Y + 2*v.X - v.Z - 1)) / 8,
		3:  -((v.Y - 1) * (v.Z - 1) * (v.Y - 2*v.X + v.Z + 1)) / 8,
		4:  -((v.Y - 1) * (v.Z + 1) * (v.Y + 2*v.X - v.Z + 1)) / 8,
		5:  -((v.Y + 1) * (v.Z + 1) * (v.Y - 2*v.X + v.Z - 1)) / 8,
		6:  ((v.Y + 1) * (v.Z + 1) * (v.Y + 2*v.X + v.Z - 1)) / 8,
		7:  ((v.Y - 1) * (v.Z + 1) * (v.Y - 2*v.X - v.Z + 1)) / 8,
		8:  -((y2 - 1) * (v.Z - 1)) / 4,
		9:  (v.X * (v.Y + 1) * (v.Z - 1)) / 2,
		10: ((y2 - 1) * (v.Z - 1)) / 4,
		11: -(v.X * (v.Y - 1) * (v.Z - 1)) / 2,
		12: ((y2 - 1) * (v.Z + 1)) / 4,
		13: -(v.X * (v.Y + 1) * (v.Z + 1)) / 2,
		14: -((y2 - 1) * (v.Z + 1)) / 4,
		15: (v.X * (v.Y - 1) * (v.Z + 1)) / 2,
		16: -((z2 - 1) * (v.Y - 1)) / 4,
		17: ((z2 - 1) * (v.Y + 1)) / 4,
		18: -((z2 - 1) * (v.Y + 1)) / 4,
		19: ((z2 - 1) * (v.Y - 1)) / 4,
		// Differentiated w.r.t Y.
		20: ((v.X - 1) * (v.Z - 1) * (2*v.Y + v.X + v.Z + 1)) / 8,
		21: -((v.X - 1) * (v.Z - 1) * (v.X - 2*v.Y + v.Z + 1)) / 8,
		22: -((v.X + 1) * (v.Z - 1) * (2*v.Y + v.X - v.Z - 1)) / 8,
		23: -((v.X + 1) * (v.Z - 1) * (2*v.Y - v.X + v.Z + 1)) / 8,
		24: -((v.X - 1) * (v.Z + 1) * (2*v.Y + v.X - v.Z + 1)) / 8,
		25: -((v.X - 1) * (v.Z + 1) * (2*v.Y - v.X + v.Z - 1)) / 8,
		26: ((v.X + 1) * (v.Z + 1) * (2*v.Y + v.X + v.Z - 1)) / 8,
		27: ((v.X + 1) * (v.Z + 1) * (2*v.Y - v.X - v.Z + 1)) / 8,
		28: -(v.Y * (v.X - 1) * (v.Z - 1)) / 2,
		29: ((x2 - 1) * (v.Z - 1)) / 4,
		30: (v.Y * (v.X + 1) * (v.Z - 1)) / 2,
		31: -((x2 - 1) * (v.Z - 1)) / 4,
		32: (v.Y * (v.X - 1) * (v.Z + 1)) / 2,
		33: -((x2 - 1) * (v.Z + 1)) / 4,
		34: -(v.Y * (v.X + 1) * (v.Z + 1)) / 2,
		35: ((x2 - 1) * (v.Z + 1)) / 4,
		36: -((z2 - 1) * (v.X - 1)) / 4,
		37: ((z2 - 1) * (v.X - 1)) / 4,
		38: -((z2 - 1) * (v.X + 1)) / 4,
		39: ((z2 - 1) * (v.X + 1)) / 4,
		// Differentiated w.r.t Z.
		40: ((v.Y - 1) * (v.X - 1) * (v.Y + v.X + 2*v.Z + 1)) / 8,
		41: -((v.Y + 1) * (v.X - 1) * (v.X - v.Y + 2*v.Z + 1)) / 8,
		42: -((v.Y + 1) * (v.X + 1) * (v.Y + v.X - 2*v.Z - 1)) / 8,
		43: -((v.Y - 1) * (v.X + 1) * (v.Y - v.X + 2*v.Z + 1)) / 8,
		44: -((v.Y - 1) * (v.X - 1) * (v.Y + v.X - 2*v.Z + 1)) / 8,
		45: -((v.Y + 1) * (v.X - 1) * (v.Y - v.X + 2*v.Z - 1)) / 8,
		46: ((v.Y + 1) * (v.X + 1) * (v.Y + v.X + 2*v.Z - 1)) / 8,
		47: ((v.Y - 1) * (v.X + 1) * (v.Y - v.X - 2*v.Z + 1)) / 8,
		48: -((y2 - 1) * (v.X - 1)) / 4,
		49: ((x2 - 1) * (v.Y + 1)) / 4,
		50: ((y2 - 1) * (v.X + 1)) / 4,
		51: -((x2 - 1) * (v.Y - 1)) / 4,
		52: ((y2 - 1) * (v.X - 1)) / 4,
		53: -((x2 - 1) * (v.Y + 1)) / 4,
		54: -((y2 - 1) * (v.X + 1)) / 4,
		55: ((x2 - 1) * (v.Y - 1)) / 4,
		56: -(v.Z * (v.Y - 1) * (v.X - 1)) / 2,
		57: (v.Z * (v.Y + 1) * (v.X - 1)) / 2,
		58: -(v.Z * (v.Y + 1) * (v.X + 1)) / 2,
		59: (v.Z * (v.Y - 1) * (v.X + 1)) / 2,
	}
}

// Quadrature returns
func (Hexa20) Quadrature() (positions []r3.Vec, weights []float64) {
	positions, weights, _ = uniformGaussQuad(3, 3, 3)
	return positions, weights
}

// String returns string representation of element type.
func (Hexa20) String() string { return "HEXA20" }

func (Hexa20) volume() float64 { return 8 }
