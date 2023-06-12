package fem

import (
	"fmt"

	"github.com/soypat/lap"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r2"
)

// AddIsoparametric3 adds isoparametric elements to the model's solid stiffness matrix.
// It calls getElement to get the element nodes and coordinates Nelem times, each with
// an incrementing element index i.
// TODO: implement arbitrary orientation of solid properties for each isoparametric element.
func (ga *GeneralAssembler) AddIsoparametric2(elemT Isoparametric2, c Constituter2D, Nelem int, getElement func(i int) (elem []int, xC r2.Vec)) error {
	const dims, dimC = 2, 3
	if elemT == nil || c == nil || getElement == nil {
		panic("nil argument to AddIsoparametric2") // This is very likely programmer error.
	}
	dofMapping, err := ga.DofMapping(elemT)
	if err != nil {
		return err
	}
	// Number of dofs per element node.
	var NelemDofsPerNode = elemT.Dofs().Count()
	if NelemDofsPerNode != dims {
		return fmt.Errorf("AddIsoparametric2: expected element to have 2 dofs per node, got %d", NelemDofsPerNode)
	}
	var (
		// Number of nodes per element.
		NnodperElem = elemT.LenNodes()
		// Number of dofs per element.
		NdofperElem = NnodperElem * NelemDofsPerNode
		// Element stiffness matrix.
		Ke = mat.NewDense(NdofperElem, NdofperElem, nil)
		// number of columns in Compliance x NdofPerNode*nodesperelement
		B = mat.NewDense(dimC, NdofperElem, nil)
		// Differentiated form functions
		dNxy = mat.NewDense(dims, NnodperElem, nil)
		// Number of dofs per node in model.
		NmodelDofsPerNode = ga.dofs.Count()
		// Quadrature integration points.
		upg, wpg = elemT.Quadrature()
	)
	C, err := c.Constitutive2D()
	if err != nil {
		return err
	}
	if len(upg) == 0 || len(upg) != len(wpg) {
		return fmt.Errorf("bad quadrature result from isoparametric element")
	}
	if r, c := C.Dims(); r != dimC || c != dimC {
		return fmt.Errorf("expected constitutive matrix to be 3x3, got %dx%d", r, c)
	}
	Cd := mat.NewDense(dimC, dimC, nil)
	Cd.Copy(C)

	// Calculate form functions evaluated at integration points.
	N := make([]*mat.VecDense, len(upg))
	dN := make([]*mat.Dense, len(upg))
	for ipg, pg := range upg {
		N[ipg] = mat.NewVecDense(NnodperElem, elemT.Basis(pg))
		dN[ipg] = mat.NewDense(dims, NnodperElem, elemT.BasisDiff(pg))
	}
	jac := mat.NewDense(dims, dims, nil)
	elemNodBacking := make([]float64, dims*NnodperElem)
	elemDofs := make([]int, NmodelDofsPerNode*NnodperElem)
	elemNod := mat.NewDense(NnodperElem, dims, elemNodBacking)
	aux1 := mat.NewDense(NdofperElem, dimC, nil)
	aux2 := mat.NewDense(NdofperElem, NdofperElem, nil)
	NvalPerElem := NdofperElem * NdofperElem
	I, J := make([]int, NvalPerElem*Nelem), make([]int, NvalPerElem*Nelem)
	V := make([]float64, NvalPerElem*Nelem)

	for iele := 0; iele < Nelem; iele++ {
		Ke.Zero()
		element, x := getElement(iele)
		if len(element) != NnodperElem {
			return fmt.Errorf("element #%d of %d nodes expected to be of %d nodes", iele, len(element), NnodperElem)
		}
		if x != (r2.Vec{}) {
			return fmt.Errorf("arbitrary constitutive orientation not implemented yet")
		}
		storeElemNode2d(elemNodBacking, ga.nodes, element)
		storeElemDofs(elemDofs, element, dofMapping, NmodelDofsPerNode)
		for ipg := range upg {
			dNi := dN[ipg]
			jac.Mul(dNi, elemNod)
			dJac := mat.Det(jac)
			if dJac < 0 {
				// return fmt.Errorf("negative determinant of jacobian of element #%d, Check node ordering", iele)
			} else if dJac < 1e-12 {
				return fmt.Errorf("zero determinant of jacobian of element #%d, Check element shape for bad aspect ratio", iele)
			}
			err := dNxy.Solve(jac, dNi)
			if err != nil {
				return fmt.Errorf("error calculating element #%d form factor: %s", iele, err)
			}

			for i := 0; i < NnodperElem; i++ {
				dNxy0 := dNxy.At(0, i)
				dNxy1 := dNxy.At(1, i)
				B.Set(0, i*dims, dNxy0)
				B.Set(1, i*dims+1, dNxy1)
				B.Set(2, i*dims, dNxy1)
				B.Set(2, i*dims+1, dNxy0)
			}
			// Ke = Ke + Bᵀ*C*B * weight*det(J)
			aux1.Mul(B.T(), Cd)
			aux2.Mul(aux1, B)
			aux2.Scale(dJac*wpg[ipg], aux2)
			Ke.Add(Ke, aux2)
		}
		offset := iele * NvalPerElem
		assembleElement(V[offset:], I[offset:], J[offset:], elemDofs, Ke)
	}
	ga.ksolid.Accumulate(lap.SparseAccum{I: I, J: J, V: V})
	return nil
}

// laptmat is a lap.Matrix wrapper that implements Gonum's Matrix interface.
type lapmat struct {
	lap.Matrix
}

func (l lapmat) T() mat.Matrix {
	if ter, ok := l.Matrix.(mat.Transpose); ok {
		return ter.T()
	}
	return lapmat{
		Matrix: lap.T(l.Matrix),
	}
}
