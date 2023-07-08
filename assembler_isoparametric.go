package fem

import (
	"errors"
	"fmt"
	"math"

	"github.com/soypat/lap"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

type quadratureCallback = func(elemIdx int, elemNodes []float64, elemDofs []int) error

// IntegrateIsoparametric is a general purpose for-each isoparametric element iterator.
func (ga *GeneralAssembler) ForEachElement(elemT Element, spatialDimsPerNode, Nelem int, getElement func(i int) []int, quadCB quadratureCallback) error {
	if elemT == nil || getElement == nil || quadCB == nil {
		panic("nil argument to AddIsoparametric2") // This is very likely programmer error.
	} else if spatialDimsPerNode < 1 || spatialDimsPerNode > 3 {
		panic("dimsPerNode must be 1, 2 or 3")
	}
	dofMapping, err := ga.DofMapping(elemT)
	if err != nil {
		return err
	}
	var (
		// Number of nodes per element.
		NnodperElem = elemT.LenNodes()
		// Number of dofs per node in model.
		NmodelDofsPerNode = ga.dofs.Count()
	)

	elemNodBacking := make([]float64, spatialDimsPerNode*NnodperElem)
	elemDofs := make([]int, NmodelDofsPerNode*NnodperElem)
	for iele := 0; iele < Nelem; iele++ {
		element := getElement(iele)
		if len(element) != NnodperElem {
			return fmt.Errorf("element #%d of %d nodes expected to be of %d nodes", iele, len(element), NnodperElem)
		}
		storeElemNode(elemNodBacking, ga.nodes, element, spatialDimsPerNode)
		storeElemDofs(elemDofs, element, dofMapping, NmodelDofsPerNode)
		err := quadCB(iele, elemNodBacking, elemDofs)
		if err != nil {
			return err
		}
	}
	return nil
}

// AddIsoparametric adds isoparametric elements to the model's solid stiffness matrix.
// It calls getElement to get the element nodes and coordinates Nelem times, each with
// an incrementing element index i.
//
// Arbitrary orientation of solid properties for each isoparametric element is not yet implemented.
func (ga *GeneralAssembler) AddIsoparametric(elemT Isoparametric, c IsoConstituter, Nelem int, getElement func(i int) (elem []int, xC, yC r3.Vec)) error {
	C, err := c.Constitutive()
	if err != nil {
		return err
	}
	dimC, _c := C.Dims()
	if _c != dimC {
		return fmt.Errorf("expected constitutive matrix to be square, got %dx%d", dimC, c)
	}

	var (
		// Number of integration dimensions per node. Usually spatial, so 2 for 2D problems and 3 for 3D problems.
		NdimsPerNode = len(elemT.BasisDiff(r3.Vec{})) / elemT.LenNodes()
		// Number of dofs per node. These contain the field variables.
		// For example, for a 2D displacement problem these are the x and y displacements, so equal to 2.
		// For a thermal problem there is always only 1 dof per node for the temperature, regardless of the number of spatial dimensions.
		NdofsPerNode = elemT.Dofs().Count() //
		// Number of nodes per element.
		NnodperElem = elemT.LenNodes()
		// Number of dofs per element.
		NdofperElem = NnodperElem * NdofsPerNode
		// Element stiffness matrix. The results of integrating the element's stiffness matrix over the element's domain.
		Ke = mat.NewDense(NdofperElem, NdofperElem, nil)
		// number of columns in Compliance x NdofPerNode*nodesperelement
		B = mat.NewDense(dimC, NdofperElem, nil)
		// Differentiated form functions with respect to the integration coordinates.
		dNxy = mat.NewDense(NdimsPerNode, NnodperElem, nil)
		// Quadrature integration points.
		upg, wpg = elemT.Quadrature()
		// Compliance matrix as a Dense matrix.
		Cd = mat.NewDense(dimC, dimC, nil)
	)
	if len(upg) == 0 || len(upg) != len(wpg) {
		return fmt.Errorf("bad quadrature result from isoparametric element")
	}
	Cd.Copy(C)

	// Calculate form functions evaluated at integration points.
	Npg := make([]*mat.VecDense, len(upg))
	dNpg := make([]*mat.Dense, len(upg))
	for ipg, pg := range upg {
		Npg[ipg] = mat.NewVecDense(NnodperElem, elemT.Basis(pg))
		dNpg[ipg] = mat.NewDense(NdimsPerNode, NnodperElem, elemT.BasisDiff(pg))
	}
	// Allocate memory for auxiliary matrices.
	jac := mat.NewDense(NdimsPerNode, NdimsPerNode, nil)

	aux1 := mat.NewDense(NdofperElem, dimC, nil)
	aux2 := mat.NewDense(NdofperElem, NdofperElem, nil)

	NvalPerElem := NdofperElem * NdofperElem
	spac := lap.NewSparseAccum(NvalPerElem * Nelem)
	var x, y r3.Vec
	subGetElement := func(i int) (elem []int) {
		elem, x, y = getElement(i)
		return elem
	}
	err = ga.ForEachElement(elemT, NdimsPerNode, Nelem, subGetElement, func(iele int, elemNodBacking []float64, elemDofs []int) error {
		if x != (r3.Vec{}) || y != (r3.Vec{}) {
			return errors.New("arbitrary constitutive orientation not implemented yet")
		}
		Ke.Zero()
		elemNod := mat.NewDense(NnodperElem, NdimsPerNode, elemNodBacking)
		for ipg := range upg {
			dN := dNpg[ipg]
			jac.Mul(dN, elemNod)
			dJac := mat.Det(jac)
			if dJac < 0 {
				return fmt.Errorf("negative determinant of jacobian of element #%d, Check node ordering", iele)
			} else if dJac < 1e-12 {
				return fmt.Errorf("zero determinant of jacobian of element #%d, Check element shape for bad aspect ratio", iele)
			}
			err := dNxy.Solve(jac, dN)
			if err != nil {
				return fmt.Errorf("error calculating element #%d form factor: %s", iele, err)
			}
			scale := c.SetStrainDisplacementMatrix(B, elemNod, dNxy, Npg[ipg])
			if math.IsNaN(scale) {
				return fmt.Errorf("NaN scale value returned by SetStrainDisplacementMatrix at element #%d, quad %d", iele, ipg)
			}
			// Ke = Ke + Bᵀ*C*B * weight*det(J)
			aux1.Mul(B.T(), Cd)
			aux2.Mul(aux1, B)
			aux2.Scale(dJac*wpg[ipg]*scale, aux2)
			Ke.Add(Ke, aux2)
		}
		offset := iele * NvalPerElem
		assembleElement(spac.V[offset:], spac.I[offset:], spac.J[offset:], elemDofs, Ke)
		return nil
	})
	if err != nil {
		return err
	}
	ga.ksolid.Accumulate(spac)
	return nil
}
