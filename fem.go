package fem

import (
	"strconv"

	"gonum.org/v1/gonum/spatial/r3"
)

// Isoparametric3 is a 3D element that employs the isoparametric formulation.
type Isoparametric3 interface {
	// LenNodes returns the number of nodes in the element.
	LenNodes() int

	// BasisNodes returns the positions of the nodes relative to the origin
	// of the element. Returned slice is of length LenNodes.
	IsoparametricNodes() []r3.Vec

	// Quadrature returns the desired quadrature integration positions and
	// respective weights.
	Quadrature() (positions []r3.Vec, weights []float64)

	// Basis returns the form functions of the element evaluated at a position
	// relative to the origin of the element. Returned slice is of length LenNodes.
	Basis(r3.Vec) []float64

	// BasisDiff returns the differentiated form functions of the element
	// evaluated at a position relative to the origin of the element.
	// The result first contains the form functions differentiated
	// with respect to X, then Y and finally Z.
	//      [ dN/dx ]
	//  N = | dN/dy |   (row major form)
	//      [ dN/dz ]
	// Suggested way of initializing matrix:
	//  dN := mat.NewDense(3, e.LenNodes(), e.BasisDiff(v))
	// Returned slice is of length LenNodes*3.
	BasisDiff(r3.Vec) []float64

	// Dofs returns the degrees of freedom corresponding to each of the element's nodes.
	Dofs() DofsFlag
}

type DofsFlag int

const (
	// DofPosX corresponds to X position degree of freedom.
	DofPosX DofsFlag = 1 << iota
	// DofPosY corresponds to Y position degree of freedom.
	DofPosY
	// DofPosZ corresponds to Z position degree of freedom.
	DofPosZ
	// DofRotX corresponds to rotational degree of freedom about X axis.
	DofRotX
	// DofRotY corresponds to rotational degree of freedom about Y axis.
	DofRotY
	// DofRotZ corresponds to rotational degree of freedom about Z axis.
	DofRotZ
	maxDofsPerNode = iota
)

// Common degree of freedom flag definitions.
const (
	DofU = DofPosX | DofPosY | DofPosZ
)

// Count returns the number of dofs set in d.
func (d DofsFlag) Count() int {
	c := 0
	for i := 0; i < maxDofsPerNode; i++ {
		c += int(d>>i) & 1
	}
	return c
}

// Has returns true if d has all of q's dofs set.
func (d DofsFlag) Has(q DofsFlag) bool {
	return d&q == q
}

// String returns a human readable representation of which dofs are set in d.
func (d DofsFlag) String() (s string) {
	for i := 0; i < maxDofsPerNode; i++ {
		if d.Has(1 << i) {
			s += strconv.Itoa(i)
		} else {
			s += "-"
		}
	}
	return s
}
