package fem

import (
	"errors"
	"fmt"
	"math/bits"
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

// Isoparametric is a 3D/2D/1D element that employs the isoparametric formulation.
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
	return bits.OnesCount16(uint16(d))
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

func forEachElement(modelDofs DofsFlag, nodes []r3.Vec, elemT Element, spatialDimsPerNode, Nelem int, getElement func(i int) []int, elemCallback elementDofCallback) error {
	if elemT == nil || getElement == nil || elemCallback == nil {
		panic("nil argument to AddIsoparametric2") // This is very likely programmer error.
	} else if spatialDimsPerNode < 1 || spatialDimsPerNode > 3 {
		return errors.New("dimsPerNode must be 1, 2 or 3")
	}

	dofMapping, err := mapdofs(modelDofs, elemT.Dofs())
	if err != nil {
		return err
	}
	var (
		// Number of nodes per element.
		NnodperElem = elemT.LenNodes()
		// Number of dofs per node in model.
		NmodelDofsPerNode = modelDofs.Count()
	)

	elemNodBacking := make([]float64, spatialDimsPerNode*NnodperElem)
	elemDofs := make([]int, NmodelDofsPerNode*NnodperElem)
	for iele := 0; iele < Nelem; iele++ {
		element := getElement(iele)
		if len(element) != NnodperElem {
			return fmt.Errorf("element #%d of %d nodes expected to be of %d nodes", iele, len(element), NnodperElem)
		}
		storeElemNode(elemNodBacking, nodes, element, spatialDimsPerNode)
		storeElemDofs(elemDofs, element, dofMapping, NmodelDofsPerNode)
		err := elemCallback(iele, elemNodBacking, elemDofs)
		if err != nil {
			return err
		}
	}
	return nil
}

func mapdofs(modelDofs, elemDofs DofsFlag) (dofs []int, err error) {
	numModelDofs := modelDofs.Count()
	if !modelDofs.Has(elemDofs) {
		return nil, fmt.Errorf("model dofs %s does not contain all element dofs %s", modelDofs.String(), elemDofs.String())
	}
	dofMapping := make([]int, elemDofs.Count())
	idm := 0
	for i := 0; i < numModelDofs; i++ {
		if modelDofs.Has(1<<i) && elemDofs.Has(1<<i) {
			dofMapping[idm] = i
			idm++
		}
	}
	if idm == 0 || len(dofMapping) == 0 {
		return nil, errors.New("element has empty dof mapping")
	}
	return dofMapping, nil
}

// func ShapeFunQuadratures[V expmat.Vector, M expmat.Matrix](e Isoparametric) ([]V, []M) {
// 	NdimsPerNode := len(e.BasisDiff(r3.Vec{})) / e.LenNodes()
// 	Nnodes := e.LenNodes()
// 	upg, _ := e.Quadrature()
// 	shape := make([]V, len(upg))
// 	shapeDiff := make([]M, len(upg))
// 	for ipg, pos := range upg {
// 		shapeDiff[ipg] = *new(M)
// 		shape[ipg] = *new(V)
// 		shape[ipg].ReuseAsVec(Nnodes)
// 		shapeDiff[ipg].ReuseAs(NdimsPerNode, Nnodes)
// 		N := e.Basis(pos)
// 		dN := e.BasisDiff(pos)
// 		for n := 0; n < Nnodes; n++ {
// 			shape[ipg].SetVec(n, N[n])
// 			for i := 0; i < NdimsPerNode; i++ {
// 				shapeDiff[ipg].Set(i, n, dN[n*NdimsPerNode+i])
// 			}
// 		}
// 	}
// 	return shape, nil
// }
