package elements

import (
	"strconv"

	"github.com/soypat/go-fem"
	"gonum.org/v1/gonum/spatial/r2"
)

// Quad8 is the 8 node 2D quadrilateral element, also known as Serendipity element
// due to the serendipitous nature of its integration points which make it unlike other 2D elements.
type Quad8 struct {
	// QuadratureOrder is the order (a.k.a degree) of the quadrature used for integration.
	// If zero a default value of 3 is used (3x3 gauss quadrature).
	QuadratureOrder int
}

var _ fem.Isoparametric2 = Quad8{}

// Dofs returns the degrees of freedom of the nodes of the element.
func (Quad8) Dofs() fem.DofsFlag {
	return fem.DofPosX | fem.DofPosY
}

// LenNodes returns the number of nodes of the element.
func (Quad8) LenNodes() int { return 8 }

// IsoparametricNodes returns the positions of the nodes relative to the origin of the element.
func (Quad8) IsoparametricNodes() []r2.Vec {
	return []r2.Vec{
		0: {X: -1, Y: -1},
		1: {X: 1, Y: -1},
		2: {X: 1, Y: 1},
		3: {X: -1, Y: 1},
		4: {X: 0, Y: -1},
		5: {X: 1, Y: 0},
		6: {X: 0, Y: 1},
		7: {X: -1, Y: 0},
	}
}

// Basis returns the form functions of the Quad8 element evaluated at v.
func (Quad8) Basis(v r2.Vec) []float64 {
	x, y := v.X, v.Y
	return []float64{
		0.25 * (1 - x) * (1 - y) * (-1 - x - y),
		0.25 * (1 + x) * (1 - y) * (-1 + x - y),
		0.25 * (1 + x) * (1 + y) * (-1 + x + y),
		0.25 * (1 - x) * (1 + y) * (-1 - x + y),
		0.5 * (1 - x) * (1 + x) * (1 - y),
		0.5 * (1 + x) * (1 + y) * (1 - y),
		0.5 * (1 - x) * (1 + x) * (1 + y),
		0.5 * (1 - x) * (1 - y) * (1 + y),
	}
}

// BasisDiff returns the differentiated form functions of the Quad8 element evaluated at v.
func (Quad8) BasisDiff(v r2.Vec) []float64 {
	x, y := v.X, v.Y
	return []float64{
		// w.r.t X
		-(x/4-1./4)*(y-1) - ((y-1)*(x+y+1))/4,
		((y-1)*(y-x+1))/4 - (x/4+1./4)*(y-1),
		(x/4+1./4)*(y+1) + ((y+1)*(x+y-1))/4,
		(x/4-1./4)*(y+1) + ((y+1)*(x-y+1))/4,
		((x+1)*(y-1))/2 + (x/2-1./2)*(y-1),
		-((y - 1) * (y + 1)) / 2,
		-((x+1)*(y+1))/2 - (x/2-1./2)*(y+1),
		((y - 1) * (y + 1)) / 2,
		// w.r.t Y
		-(x/4-1./4)*(y-1) - (x/4-1./4)*(x+y+1),
		(x/4+1./4)*(y-x+1) + (x/4+1./4)*(y-1),
		(x/4+1./4)*(y+1) + (x/4+1./4)*(x+y-1),
		(x/4-1./4)*(x-y+1) - (x/4-1./4)*(y+1),
		(x/2 - 1./2) * (x + 1),
		-(x/2+1./2)*(y-1) - (x/2+1./2)*(y+1),
		-(x/2 - 1./2) * (x + 1),
		(x/2-1./2)*(y-1) + (x/2-1./2)*(y+1),
	}
}

// Quadrature returns the quadrature nodes and weights of the element.
func (q8 Quad8) Quadrature() ([]r2.Vec, []float64) {
	quad := q8.order()
	pos, w, err := uniformGaussQuad2d(quad, quad)
	if err != nil {
		panic(err)
	}
	return pos, w
}

func (q8 Quad8) order() int {
	quad := q8.QuadratureOrder
	if quad <= 0 {
		quad = 3
	}
	return quad
}

// String returns a string representation of the element.
func (q8 Quad8) String() string { return "QUAD8(order=" + strconv.Itoa(q8.order()) + ")" }

func (Quad8) area() float64 { return 4 }
