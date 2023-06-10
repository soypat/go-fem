package elements

import (
	"strconv"

	"github.com/soypat/go-fem"
	"gonum.org/v1/gonum/spatial/r2"
)

// Quad4 is the 4 node 2D quadrilateral element.
type Quad4 struct {
	// QuadratureOrder is the order (a.k.a degree) of the quadrature used for integration.
	// If zero a default value of 2 is used (2x2 gauss quadrature).
	QuadratureOrder int
}

var _ fem.Isoparametric2 = Quad4{}

// Dofs returns the degrees of freedom of the nodes of the element.
func (Quad4) Dofs() fem.DofsFlag {
	return fem.DofPosX | fem.DofPosY
}

// LenNodes returns the number of nodes of the element.
func (Quad4) LenNodes() int { return 4 }

// IsoparametricNodes returns the positions of the nodes relative to the origin of the element.
func (Quad4) IsoparametricNodes() []r2.Vec {
	return []r2.Vec{
		0: {X: -1, Y: -1},
		1: {X: 1, Y: -1},
		2: {X: 1, Y: 1},
		3: {X: -1, Y: 1},
	}
}

// Basis returns the form functions of the Quad4 element evaluated at v.
func (Quad4) Basis(v r2.Vec) []float64 {
	x, y := v.X, v.Y
	return []float64{
		0.25 * (1 - x) * (1 - y),
		0.25 * (1 + x) * (1 - y),
		0.25 * (1 + x) * (1 + y),
		0.25 * (1 - x) * (1 + y),
	}
}

// BasisDiff returns the differentiated form functions of the Quad4 element evaluated at v.
func (Quad4) BasisDiff(v r2.Vec) []float64 {
	x, y := v.X, v.Y
	return []float64{
		// w.r.t X
		y/4 - 1./4,
		1./4 - y/4,
		y/4 + 1./4,
		-y/4 - 1./4,
		// w.r.t Y
		x/4 - 1./4,
		-x/4 - 1./4,
		x/4 + 1./4,
		1./4 - x/4,
	}
}

// Quadrature returns the quadrature nodes and weights of the element.
func (q8 Quad4) Quadrature() ([]r2.Vec, []float64) {
	quad := q8.order()
	pos, w, err := uniformGaussQuad2d(quad, quad)
	if err != nil {
		panic(err)
	}
	return pos, w
}

func (q4 Quad4) order() int {
	quad := q4.QuadratureOrder
	if quad <= 0 {
		quad = 2
	}
	return quad
}

// String returns a string representation of the element.
func (q4 Quad4) String() string { return "QUAD4(order=" + strconv.Itoa(q4.order()) + ")" }

func (Quad4) area() float64 { return 4 }
