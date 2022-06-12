package elements

import (
	"github.com/soypat/go-fem"

	"gonum.org/v1/gonum/spatial/r3"
)

// Tetra4 is the 3D linear strain tetrahedral element of 4 nodes.
type Tetra4 struct{}

var _ fem.Isoparametric3 = Tetra4{}

func (Tetra4) Dofs() fem.DofsFlag {
	return fem.DofU
}

func (Tetra4) LenNodes() int { return 4 }

// IsoparametricNodes returns the positions of the nodes relative to the origin of the element.
func (Tetra4) IsoparametricNodes() []r3.Vec {
	return []r3.Vec{
		0: {X: 0, Y: 0, Z: 0},
		1: {X: 1, Y: 0, Z: 0},
		2: {X: 0, Y: 1, Z: 0},
		3: {X: 0, Y: 0, Z: 1},
	}
}

// Basis returns the form functions of the Tetra10 element evaluated at v.
func (Tetra4) Basis(v r3.Vec) []float64 {
	return []float64{1 - v.X - v.Y - v.Z, v.X, v.Y, v.Z}

}

func (Tetra4) BasisDiff(v r3.Vec) []float64 {
	return []float64{
		-1, 1, 0, 0, // w.r.t X
		-1, 0, 1, 0, // w.r.t Y
		-1, 0, 0, 1, // w.r.t Z
	}
}

func (Tetra4) Quadrature() (positions []r3.Vec, weights []float64) {
	return []r3.Vec{{X: .25, Y: .25, Z: .25}}, []float64{1}
}
