package fem

import (
	"errors"
	"fmt"

	"github.com/soypat/go-fem/internal/expmat"
	"github.com/soypat/lap"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

// GeneralAssembler is a general purpose stiffness matrix assembler.
type GeneralAssembler struct {
	// Stiffness matrix of modelled solid.
	ksolid lap.Sparse
	nodes  []r3.Vec
	dofs   DofsFlag
}

// NewSymAssembler initializes a GeneralAssembler ready for use.
func NewGeneralAssembler(nodes []r3.Vec, modelDofs DofsFlag) *GeneralAssembler {
	totalDofs := len(nodes) * modelDofs.Count()
	return &GeneralAssembler{
		ksolid: *lap.NewSparse(totalDofs, totalDofs),
		dofs:   modelDofs,
		nodes:  nodes,
	}
}

// Ksolid returns the stiffness matrix of the solid.
func (ga *GeneralAssembler) Ksolid() *lap.Sparse { return &ga.ksolid }

// assembleElement stores the element stiffness matrix Ke into V data for the corresponding
// I and J indices to the global stiffness matrix.
func assembleElement(V []float64, I, J, elemDofs []int, Ke *mat.Dense) {
	_, c := Ke.Dims()
	for i, ei := range elemDofs {
		ic := i * c
		row := Ke.RawRowView(i)
		copy(V[ic:], row)
		copy(J[ic:], elemDofs)
		// TODO(soypat): could this be a copy or a optimized iteration?
		for j := range elemDofs {
			ix := ic + j
			I[ix] = ei
		}
	}
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
	if idm == 0 || len(dofMapping) == 0 {
		return nil, errors.New("element has empty dof mapping")
	}
	return dofMapping, nil
}

func storeElemNode(dst []float64, allNodes []r3.Vec, elem []int, dims int) {
	if len(dst) != dims*len(elem) {
		panic("bad length")
	}
	switch dims {
	case 1:
		for i := range elem {
			dst[i] = allNodes[elem[i]].X
		}
	case 2:
		for i := range elem {
			offset := i * 2
			n := allNodes[elem[i]]
			dst[offset] = n.X
			dst[offset+1] = n.Y
		}
	case 3:
		for i := range elem {
			offset := i * 3
			n := allNodes[elem[i]]
			dst[offset] = n.X
			dst[offset+1] = n.Y
			dst[offset+2] = n.Z
		}
	default:
		panic("unhandled dimension")
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
		if x != (r3.Vec{}) || y != (r3.Vec{}) {
			return errors.New("element orientation not supported yet")
		}
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

		if x != (r3.Vec{}) && y != (r3.Vec{}) {
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
		}
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
	if i > b.rep*r || j > b.rep*c {
		panic("bad access")
	}
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
