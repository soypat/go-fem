package elements

import (
	"strconv"

	"github.com/soypat/go-fem"
	"gonum.org/v1/gonum/spatial/r2"
)

// Triangle3 is the 3 node 2D triangle element.
type Triangle3 struct {
	// QuadratureOrder is the degree of the quadrature used for integration.
	// If zero a default value of 1 is used (1 node gauss quadrature).
	QuadratureOrder int
}

var _ fem.Isoparametric2 = Triangle3{}

// Dofs returns the degrees of freedom of the nodes of the element.
func (Triangle3) Dofs() fem.DofsFlag {
	return fem.DofPosX | fem.DofPosY
}

// LenNodes returns the number of nodes of the element.
func (Triangle3) LenNodes() int { return 3 }

// IsoparametricNodes returns the positions of the nodes relative to the origin of the element.
func (Triangle3) IsoparametricNodes() []r2.Vec {
	const b, h = 1, 1
	return []r2.Vec{
		0: {X: 0, Y: 0},
		1: {X: b, Y: 0},
		2: {X: 0, Y: h},
	}
}

// Basis returns the form functions of the Triangle3 element evaluated at v.
func (Triangle3) Basis(v r2.Vec) []float64 {
	x, y := v.X, v.Y
	return []float64{
		1 - x - y,
		x,
		y,
	}
}

// BasisDiff returns the differentiated form functions of the Triangle3 element evaluated at v.
func (Triangle3) BasisDiff(v r2.Vec) []float64 {
	return []float64{
		-1,
		1,
		0,
		-1,
		0,
		1,
	}
}

// Quadrature returns the quadrature integration nodes and weights of the element.
func (t3 Triangle3) Quadrature() (nodes []r2.Vec, weights []float64) {
	return getTriangleQuads(t3.order())
}

func (t3 Triangle3) order() int {
	quad := t3.QuadratureOrder
	if quad <= 0 {
		quad = 1
	}
	return quad
}

func (t3 Triangle3) String() string { return "TRI6(order=" + strconv.Itoa(t3.order()) + ")" }

func (Triangle3) area() float64 { return 0.5 }

// Triangle6 is the 6 node 2D triangle element. The 3 nodes besides the corners
// are located at the middle of the edges.
type Triangle6 struct {
	// QuadratureOrder is the degree of the quadrature used for integration.
	// If zero a default value of 2 is used (3 node gauss quadrature).
	QuadratureOrder int
}

var _ fem.Isoparametric2 = Triangle6{}

// Dofs returns the degrees of freedom of the nodes of the element.
func (Triangle6) Dofs() fem.DofsFlag {
	return fem.DofPosX | fem.DofPosY
}

// LenNodes returns the number of nodes of the element.
func (Triangle6) LenNodes() int { return 6 }

// IsoparametricNodes returns the positions of the nodes relative to the origin of the element.
func (Triangle6) IsoparametricNodes() []r2.Vec {
	const b, h = 1, 1
	return []r2.Vec{
		0: {X: 0, Y: 0},
		1: {X: b, Y: 0},
		2: {X: 0, Y: h},
		3: {X: b / 2, Y: h / 2},
		4: {X: 0, Y: h / 2},
		5: {X: b / 2, Y: 0},
	}
}

// Basis returns the form functions of the Triangle6 element evaluated at v.
func (Triangle6) Basis(v r2.Vec) []float64 {
	x, y := v.X, v.Y
	return []float64{
		1 - 3*x - 3*y + 2*x*x + 4*x*y + 2*y*y,
		4*x - 4*x*x - 4*x*y,
		-x + 2*x*x,
		4 * x * y,
		-y + 2*y*y,
		4*y - 4*x*y - 4*y*y,
	}
}

// BasisDiff returns the differentiated form functions of the Triangle6 element evaluated at v.
func (Triangle6) BasisDiff(v r2.Vec) []float64 {
	x, y := v.X, v.Y
	return []float64{
		-3 + 4*x + 4*y,
		4 - 8*x - 4*y,
		-1 + 4*x,
		4 * y,
		-4 * y,
		-4 + 4*x + 8*y,
		4 * x,
		-4 * x,
	}
}

// Quadrature returns the quadrature integration nodes and weights of the element.
func (t6 Triangle6) Quadrature() (nodes []r2.Vec, weights []float64) {
	return getTriangleQuads(t6.order())
}

func (t6 Triangle6) order() int {
	quad := t6.QuadratureOrder
	if quad <= 0 {
		quad = 2
	}
	return quad
}

func (t6 Triangle6) String() string { return "TRI6(order=" + strconv.Itoa(t6.order()) + ")" }

func (Triangle6) area() float64 { return 0.5 }

func getTriangleQuads(order int) (nodes []r2.Vec, weights []float64) {
	if order < 1 {
		panic("bad quadrature order; bug in go-fem")
	}
	if order > len(triangleQuads) {
		panic("triangle quadrature order not implemented")
	}
	// Clone slices.
	quads := triangleQuads[order-1]
	nodes = append([]r2.Vec{}, quads.Nodes...)
	weights = append([]float64{}, quads.Weights...)
	for i := range weights {
		weights[i] *= 0.5 // Denormalize weights.
	}
	return nodes, weights
}

// https://mathsfromnothing.au/triangle-quadrature-rules/?i=1
var triangleQuads = []struct {
	Weights []float64
	Nodes   []r2.Vec
}{
	0: { // 1 node
		Weights: []float64{1},
		Nodes:   []r2.Vec{{X: 1.0 / 3, Y: 1.0 / 3}},
	},
	1: { // Order 2
		Weights: []float64{1.0 / 3, 1.0 / 3, 1.0 / 3},
		Nodes: []r2.Vec{
			{X: 1.0 / 6, Y: 2.0 / 3},
			{X: 1.0 / 6, Y: 1.0 / 6},
			{X: 2.0 / 3, Y: 1.0 / 6},
		},
	},
	2: { // Order 3
		Weights: []float64{-0.5625, 0.520833333333333, 0.520833333333333, 0.520833333333333},
		Nodes: []r2.Vec{
			{X: 1.0 / 3, Y: 1.0 / 3},
			{X: 0.2, Y: 0.6},
			{X: 0.2, Y: 0.2},
			{X: 0.6, Y: 0.2},
		},
	},
	3: { // Order 4.
		Weights: []float64{0.223381589678011, 0.223381589678011, 0.223381589678011, 0.109951743655322, 0.109951743655322, 0.109951743655322},
		Nodes: []r2.Vec{
			{X: 0.445948490915965, Y: 0.108103018168070},
			{X: 0.445948490915965, Y: 0.445948490915965},
			{X: 0.108103018168070, Y: 0.445948490915965},
			{X: 0.091576213509771, Y: 0.816847572980459},
			{X: 0.091576213509771, Y: 0.091576213509771},
			{X: 0.816847572980459, Y: 0.091576213509771},
		},
	},
}
