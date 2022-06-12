package fem

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

// GeneralAssembler is a general purpose stiffness matrix assembler.
type GeneralAssembler struct {
	// Stiffness matrix of modelled solid.
	ksolid mat.Dense
	dofs   DofsFlag
	nodes  []r3.Vec
}

// NewSymAssembler initializes a GeneralAssembler ready for use.
func NewGeneralAssembler(nodes []r3.Vec, modelDofs DofsFlag) *GeneralAssembler {
	totalDofs := len(nodes) * modelDofs.Count()
	return &GeneralAssembler{
		ksolid: *mat.NewDense(totalDofs, totalDofs, nil),
		dofs:   modelDofs,
		nodes:  nodes,
	}
}

// Ksolid returns the stiffness matrix of the solid.
func (ga *GeneralAssembler) Ksolid() mat.Matrix { return &ga.ksolid }

// AddIsoparametric3s adds isoparametric elements to the model's solid stiffness matrix.
func (ga *GeneralAssembler) AddIsoparametric3s(elemT Isoparametric3, c Constituter, Nelem int, getElement func(i int) []int) error {
	if elemT == nil || c == nil || getElement == nil {
		panic("nil argument to AddIsoparametric3s") // This is very likely programmer error.
	}
	if !ga.dofs.Has(elemT.Dofs()) {
		return fmt.Errorf("model dofs %s does not contain all element dofs %s", ga.dofs.String(), elemT.Dofs().String())
	}
	var (
		// Number of nodes per element.
		NnodperElem = elemT.LenNodes()
		// Number of dofs per element node.
		NelemDofsPerNode = elemT.Dofs().Count()
		// Number of dofs per element.
		NdofperElem = NnodperElem * NelemDofsPerNode
		// Element stiffness matrix.
		Ke = mat.NewDense(NdofperElem, NdofperElem, nil)
		// number of columns in Compliance x NdofPerNode*nodesperelement
		B = mat.NewDense(6, NdofperElem, nil)
		// Differentiated form functions
		dNxyz             = mat.NewDense(3, NnodperElem, nil)
		C                 = c.Constitutive()
		NmodelDofsPerNode = ga.dofs.Count()
		// Quadrature integration points.
		upg, wpg = elemT.Quadrature()
	)
	if len(upg) == 0 || len(upg) != len(wpg) {
		return fmt.Errorf("bad quadrature result from isoparametric element")
	}
	if r, c := C.Dims(); r != 6 || c != 6 {
		return fmt.Errorf("expected constitutive matrix to be 6x6, got %dx%d", r, c)
	}
	// Create element to model dof mapping.
	var dofMapping []int
	for i := 0; i < ga.dofs.Count(); i++ {
		if ga.dofs.Has(1<<i) && elemT.Dofs().Has(1<<i) {
			dofMapping = append(dofMapping, i)
		}
	}
	// Calculate form functions evaluated at integration points.
	N := make([]*mat.VecDense, len(upg))
	dN := make([]*mat.Dense, len(upg))
	for ipg, pg := range upg {
		N[ipg] = mat.NewVecDense(NnodperElem, elemT.Basis(pg))
		dN[ipg] = mat.NewDense(3, NnodperElem, elemT.BasisDiff(pg))
	}
	var jac r3.Mat
	elemNodBacking := make([]float64, 3*NnodperElem)
	elemDofs := make([]int, NmodelDofsPerNode*NnodperElem)
	elemNod := mat.NewDense(NnodperElem, 3, elemNodBacking)
	aux1 := mat.NewDense(NdofperElem, 6, nil)
	aux2 := mat.NewDense(NdofperElem, NdofperElem, nil)
	for iele := 0; iele < Nelem; iele++ {
		Ke.Zero()
		element := getElement(iele)
		if len(element) != NnodperElem {
			return fmt.Errorf("element #%d of %d nodes expected to be of %d nodes", iele, len(element), NnodperElem)
		}
		storeElemNode(elemNodBacking, ga.nodes, element)
		storeElemDofs(elemDofs, element, dofMapping, NmodelDofsPerNode)
		for ipg := range upg {
			dNi := dN[ipg]
			jac.Mul(dNi, elemNod)
			err := dNxyz.Solve(&jac, dNi)
			if err != nil {
				return fmt.Errorf("error calculating element #%d form factor: %s", iele, err)
			}
			for i := 0; i < NnodperElem; i++ {
				// First three rows.
				B.Set(0, i*3, dNxyz.At(0, i))
				B.Set(1, i*3+1, dNxyz.At(1, i))
				B.Set(2, i*3+2, dNxyz.At(2, i))
				// Fourth row.
				B.Set(3, i*3, dNxyz.At(1, i))
				B.Set(3, i*3+1, dNxyz.At(0, i))
				// Fifth row.
				B.Set(4, i*3+1, dNxyz.At(2, i))
				B.Set(4, i*3+2, dNxyz.At(1, i))
				// Sixth row.
				B.Set(5, i*3, dNxyz.At(2, i))
				B.Set(5, i*3+2, dNxyz.At(0, i))
			}
			// Ke = Ke + Bᵀ*C*B * weight*det(J)
			aux1.Mul(B.T(), C)
			aux2.Mul(aux1, B)
			aux2.Scale(jac.Det()*wpg[ipg], aux2)
			Ke.Add(Ke, aux2)
		}
		// Ke is a matrix of reduced dof size.
		for i := 0; i < NdofperElem; i++ {
			ei := elemDofs[i]
			for j := 0; j < NdofperElem; j++ {
				ej := elemDofs[j]
				ga.ksolid.Set(ei, ej, ga.ksolid.At(ei, ej)+Ke.At(i, j))
			}
		}
	}
	return nil
}

func storeElemNode(dst []float64, allNodes []r3.Vec, elem []int) {
	if len(dst) != 3*len(elem) {
		panic("bad length")
	}
	for i := range elem {
		offset := i * 3
		n := allNodes[elem[i]]
		dst[offset] = n.X
		dst[offset+1] = n.Y
		dst[offset+2] = n.Z
	}
}

func storeElemDofs(dst, elem, dofmapping []int, modelDofsPerNode int) {
	if len(dofmapping)*len(elem) != len(dst) {
		panic("bad length")
	}
	for inod, node := range elem {
		offset := inod * len(dofmapping)
		dofstart := modelDofsPerNode * node
		for j, dofoffset := range dofmapping {
			dst[offset+j] = dofstart + dofoffset
		}
	}
}
