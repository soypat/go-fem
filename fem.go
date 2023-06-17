package fem

import (
	"strconv"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

type Element interface {
	// LenNodes returns the number of nodes in the element.
	LenNodes() int
	// Dofs returns the degrees of freedom corresponding to each of the element's nodes.
	Dofs() DofsFlag
}

// Isoparametric is a 3D element that employs the isoparametric formulation.
type Isoparametric interface {
	Element
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
}

type Element3 interface {
	Element
	CopyK(dst *mat.Dense, elementNodes []r3.Vec) error
	SetConstitutive(c Constituter) error
}

// Constituter represents the homogenous properties of a medium
// that can then be used to model solids or other continuous field problems.
// For solids it returns the unmodified constitutive tensor (Generalized Hookes law).
// The shear modulii should not be halved.
type Constituter interface {
	Constitutive() (mat.Matrix, error)
}

type IsoConstituter interface {
	Constituter
	SetStrainDisplacementMatrix(dstB, elemNodes, dN *mat.Dense, N *mat.VecDense) (scale float64)
}

// DofsFlag holds bitwise information of degrees of freedom.
type DofsFlag uint16

const (
	// DofPosX corresponds to X position degree of freedom (0).
	DofPosX DofsFlag = 1 << iota
	// DofPosY corresponds to Y position degree of freedom (1).
	DofPosY
	// DofPosZ corresponds to Z position degree of freedom (2).
	DofPosZ
	// DofRotX corresponds to rotational degree of freedom about X axis (3).
	DofRotX
	// DofRotY corresponds to rotational degree of freedom about Y axis (4).
	DofRotY
	// DofRotZ corresponds to rotational degree of freedom about Z axis (5).
	DofRotZ
	maxDofsPerNode = iota
)

// Common degree of freedom flag definitions.
const (
	DofPos = DofPosX | DofPosY | DofPosZ
	DofRot = DofRotX | DofRotY | DofRotZ
	Dof6   = DofPos | DofRot
)

// Count returns the number of dofs set in d.
func (d DofsFlag) Count() int {
	c := 0
	for i := 0; i < maxDofsPerNode; i++ {
		c += int(d>>i) & 1
	}
	return c
}

// Has returns true if d has all of q's dofs set. It returns false if
// q has a dof set that is not set in d.
func (d DofsFlag) Has(q DofsFlag) bool {
	return d&q == q
}

// String returns a human readable representation of which dofs are set in d.
func (d DofsFlag) String() (s string) {
	if d == 0 {
		return "none"
	}
	for i := 0; d != 0; i++ {
		if d&1 != 0 {
			s += strconv.Itoa(i + 1)
		} else {
			s += "-"
		}
		d = d >> 1
	}
	return s
}
