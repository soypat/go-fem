package fem

import (
	"fmt"

	"github.com/soypat/go-fem/internal/expmat"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

// GeneralAssembler is a general purpose stiffness matrix assembler.
type GeneralAssembler struct {
	// Stiffness matrix of modelled solid.
	ksolid expmat.Sparse
	dofs   DofsFlag
	nodes  []r3.Vec
}

// NewSymAssembler initializes a GeneralAssembler ready for use.
func NewGeneralAssembler(nodes []r3.Vec, modelDofs DofsFlag) *GeneralAssembler {
	totalDofs := len(nodes) * modelDofs.Count()
	return &GeneralAssembler{
		ksolid: *expmat.NewSparse(totalDofs, totalDofs),
		dofs:   modelDofs,
		nodes:  nodes,
	}
}

// Ksolid returns the stiffness matrix of the solid.
func (ga *GeneralAssembler) Ksolid() *expmat.Sparse { return &ga.ksolid }

// AddIsoparametric3 adds isoparametric elements to the model's solid stiffness matrix.
// TODO: implement arbitrary orientation of solid properties for each isoparametric element.
func (ga *GeneralAssembler) AddIsoparametric3(elemT Isoparametric3, c Constituter, Nelem int, getElement func(i int) (elem []int, xC, yC r3.Vec)) error {
	if elemT == nil || c == nil || getElement == nil {
		panic("nil argument to AddIsoparametric3") // This is very likely programmer error.
	}
	dofMapping, err := ga.DofMapping(elemT)
	if err != nil {
		return err
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
		element, x, y := getElement(iele)
		if len(element) != NnodperElem {
			return fmt.Errorf("element #%d of %d nodes expected to be of %d nodes", iele, len(element), NnodperElem)
		}
		if x != (r3.Vec{}) || y != (r3.Vec{}) {
			return fmt.Errorf("arbitrary constitutive orientation not implemented yet")
		}
		storeElemNode(elemNodBacking, ga.nodes, element)
		storeElemDofs(elemDofs, element, dofMapping, NmodelDofsPerNode)
		for ipg := range upg {
			dNi := dN[ipg]
			jac.Mul(dNi, elemNod)
			dJac := jac.Det()
			if dJac < 0 {
				// return fmt.Errorf("negative determinant of jacobian of element #%d, Check node ordering", iele)
			} else if dJac < 1e-12 {
				return fmt.Errorf("zero determinant of jacobian of element #%d, Check element shape for bad aspect ratio", iele)
			}
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
			aux2.Scale(dJac*wpg[ipg], aux2)
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

// DofMapping creates element to model dof mapping. If model does not contain all argument elements
// dofs then an error is returned.
// The length of the slice returned is equal to the amount of dofs per element node.
func (ga *GeneralAssembler) DofMapping(e Element) ([]int, error) {
	edofs := e.Dofs()
	if !ga.dofs.Has(edofs) {
		return nil, fmt.Errorf("model dofs %s does not contain all element dofs %s", ga.dofs.String(), e.Dofs().String())
	}
	dofMapping := make([]int, edofs.Count())
	idm := 0
	for i := 0; i < ga.dofs.Count(); i++ {
		if ga.dofs.Has(1<<i) && edofs.Has(1<<i) {
			dofMapping[idm] = i
			idm++
		}
	}
	return dofMapping, nil
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

func (ga *GeneralAssembler) AddElement3(elemT Element3, c Constituter, Nelem int, getElement func(i int) (e []int, x, y r3.Vec)) error {
	if elemT == nil || c == nil || getElement == nil {
		panic("nil argument to AddElement3") // This is very likely programmer error.
	}
	dofMapping, err := ga.DofMapping(elemT)
	if err != nil {
		return err
	}
	var (
		// Number of dofs per element node.
		Nde = len(dofMapping)
		// Number of nodes per element.
		NnodPerElem = elemT.LenNodes()
		// Number of dofs per element.
		NdofPerElem      = Nde * NnodPerElem
		NdofPerNodeModel = ga.dofs.Count()
	)
	if ga.dofs>>6 != 0 {
		panic("AddElement3 currently only handles 6 rigid body motion degrees of freedom")
	}
	var T3 r3.Mat
	var rotator mat.Matrix = &expmat.SubMat{
		Rix: dofMapping,
		Cix: dofMapping,
		M: blkDiag{
			rep: 2 * NnodPerElem,
			m:   &T3,
		},
	}

	Ke := mat.NewDense(NdofPerElem, NdofPerElem, nil)
	elemNodes := make([]r3.Vec, NnodPerElem)
	elemDofs := make([]int, NdofPerElem)
	for iele := 0; iele < Nelem; iele++ {
		Ke.Zero()
		element, x, y := getElement(iele)
		if len(element) != NnodPerElem {
			return fmt.Errorf("element #%d of %d nodes expected to be of %d nodes", iele, len(element), NnodPerElem)
		}
		for i, elnod := range element {
			elemNodes[i] = ga.nodes[elnod]
		}
		err := elemT.CopyK(Ke, elemNodes)
		if err != nil {
			return err
		}
		// Rotate element stiffness matrix to match user input orientation.
		orientX := r3.Unit(x)
		orientY := r3.Unit(y)
		orientZ := r3.Cross(orientX, orientY)
		orientY = r3.Cross(orientZ, orientX) // Should be unit vector.
		T3.Set(0, 0, orientX.X)
		T3.Set(0, 1, orientX.Y)
		T3.Set(0, 2, orientX.Z)
		T3.Set(1, 0, orientY.X)
		T3.Set(1, 1, orientY.Y)
		T3.Set(1, 2, orientY.Z)
		T3.Set(2, 0, orientZ.X)
		T3.Set(2, 1, orientZ.Y)
		T3.Set(2, 2, orientZ.Z)
		Ke.Mul(rotator.T(), Ke)
		Ke.Mul(Ke, rotator)
		storeElemDofs(elemDofs, element, dofMapping, NdofPerNodeModel)
		for i := 0; i < NdofPerElem; i++ {
			ei := elemDofs[i]
			for j := 0; j < NdofPerElem; j++ {
				ej := elemDofs[j]
				ga.ksolid.Set(ei, ej, ga.ksolid.At(ei, ej)+Ke.At(i, j))
			}
		}
	}
	return nil
}

type blkDiag struct {
	rep int
	m   mat.Matrix
}

func (b blkDiag) Dims() (int, int) {
	r, c := b.m.Dims()
	return r * b.rep, c * b.rep
}

func (b blkDiag) At(i, j int) float64 {
	r, c := b.m.Dims()
	iq := i / r
	jq := j / c
	if iq != jq {
		return 0
	}
	return b.m.At(i%r, j%c)
}

func (b blkDiag) T() mat.Matrix {
	return blkDiag{
		rep: b.rep,
		m:   b.m.T(),
	}
}
