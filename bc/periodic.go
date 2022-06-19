package bc

import (
	"errors"
	"fmt"
	"unsafe"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/exp/expmat"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

func (bf BoxFeatures) PeriodicImposedDisplacements(dst []float64, rucCase int, displacement float64) error {
	rows := 0
	imposedDisp := imposedDisplacementForRUC(rucCase, 0.1)
	dx := imposedDisp.At(0, 0)
	dy := imposedDisp.At(1, 1)
	dz := imposedDisp.At(2, 2)
	dxy := imposedDisp.At(0, 1)
	dxz := imposedDisp.At(0, 2)
	dyz := imposedDisp.At(1, 2)
	modelSize := r3.Sub(bf.bounds.Max, bf.bounds.Min)
	bfc, err := bf.countFeatures()
	if err != nil {
		return err
	}
	var imposed [3]float64
	// Surface load imposition (Lagrange).
	imposed = [3]float64{dx * modelSize.X, dxy * modelSize.X, dxz * modelSize.X}
	for i := 0; i < bfc.sx; i++ {
		copy(dst[rows*3:], imposed[:])
		rows++
	}
	imposed = [3]float64{dxy * modelSize.Y, dy * modelSize.Y, dyz * modelSize.Y}
	for i := 0; i < bfc.sy; i++ {
		copy(dst[rows*3:], imposed[:])
		rows++
	}
	imposed = [3]float64{dxz * modelSize.Z, dyz * modelSize.Z, dz * modelSize.Z}
	for i := 0; i < bfc.sz; i++ {
		copy(dst[rows*3:], imposed[:])
		rows++
	}
	// Edge load imposition (Lagrange).
	imposed = [3]float64{modelSize.X*dx - modelSize.Z*dxz, modelSize.X*dxy - modelSize.Z*dyz, modelSize.X*dxz - modelSize.Z*dz}
	for i := 0; i < bfc.ne1; i++ {
		copy(dst[rows*3:], imposed[:])
		rows++
	}
	imposed = [3]float64{modelSize.X*dx + modelSize.Z*dxz, modelSize.X*dxy + modelSize.Z*dyz, modelSize.X*dxz + modelSize.Z*dz}
	for i := 0; i < bfc.ne2; i++ {
		copy(dst[rows*3:], imposed[:])
		rows++
	}
	imposed = [3]float64{modelSize.X*dx + modelSize.Y*dxy, modelSize.X*dxy + modelSize.Y*dy, modelSize.X*dxz + modelSize.Y*dyz}
	for i := 0; i < bfc.ne3; i++ {
		copy(dst[rows*3:], imposed[:])
		rows++
	}
	imposed = [3]float64{modelSize.X*dx - modelSize.Y*dxy, modelSize.X*dxy - modelSize.Y*dy, modelSize.X*dxz - modelSize.Y*dyz}
	for i := 0; i < bfc.ne4; i++ {
		copy(dst[rows*3:], imposed[:])
		rows++
	}
	imposed = [3]float64{modelSize.Y*dxy - modelSize.Z*dxz, modelSize.Y*dy - modelSize.Z*dyz, modelSize.Y*dyz - modelSize.Z*dz}
	for i := 0; i < bfc.ne5; i++ {
		copy(dst[rows*3:], imposed[:])
		rows++
	}
	imposed = [3]float64{modelSize.Y*dxy + modelSize.Z*dxz, modelSize.Y*dy + modelSize.Z*dyz, modelSize.Y*dyz + modelSize.Z*dz}
	for i := 0; i < bfc.ne6; i++ {
		copy(dst[rows*3:], imposed[:])
		rows++
	}
	// Corner load imposition (Lagrange).
	imposed = [3]float64{modelSize.X*dx + modelSize.Y*dxy - modelSize.Z*dxz, modelSize.X*dxy + modelSize.Y*dy - modelSize.Z*dyz, modelSize.X*dxz + modelSize.Y*dyz - modelSize.Z*dz}
	copy(dst[rows*3:], imposed[:])
	rows++
	imposed = [3]float64{modelSize.X*dx + modelSize.Y*dxy + modelSize.Z*dxz, modelSize.X*dxy + modelSize.Y*dy + modelSize.Z*dyz, modelSize.X*dxz + modelSize.Y*dyz + modelSize.Z*dz}
	copy(dst[rows*3:], imposed[:])
	rows++
	imposed = [3]float64{-modelSize.X*dx + modelSize.Y*dxy + modelSize.Z*dxz, -modelSize.X*dxy + modelSize.Y*dy + modelSize.Z*dyz, -modelSize.X*dxz + modelSize.Y*dyz + modelSize.Z*dz}
	copy(dst[rows*3:], imposed[:])
	rows++
	imposed = [3]float64{modelSize.X*dx - modelSize.Y*dxy + modelSize.Z*dxz, modelSize.X*dxy - modelSize.Y*dy + modelSize.Z*dyz, modelSize.X*dxz - modelSize.Y*dyz + modelSize.Z*dz}
	copy(dst[rows*3:], imposed[:])
	return nil
}

func imposedDisplacementForRUC(rucCase int, displacement float64) *mat.Dense {
	var dx, dy, dz, dxy, dxz, dyz float64
	switch rucCase {
	case 0:
		dx = displacement
	case 1:
		dy = displacement
	case 2:
		dz = displacement
	case 3:
		dxy = displacement / 2
	case 4:
		dxz = displacement / 2
	case 5:
		dyz = displacement / 2
	default:
		panic("invalid RUC case")
	}
	return mat.NewDense(3, 3, []float64{
		dx, dxy, dxz,
		dxy, dy, dyz,
		dxz, dyz, dz,
	})
}

type BoxFeatures struct {
	nodes  []r3.Vec
	bounds r3.Box
	// RUC surfaces
	sx, sX, sy, sY, sz, sZ []int
	// RUC edges.
	exy, eXy, exY, eXY, exz, eXz, exZ, eXZ, eyz, eYz, eyZ, eYZ []int
	// RUC corners.
	cxyz, cXyz, cxYz, cXYz, cxyZ, cXyZ, cxYZ, cXYZ int
}

type boxFeatureCount struct {
	// Number of surface nodes.
	sx, sy, sz int
	// Number of edge nodes.
	ne1, ne2, ne3, ne4, ne5, ne6 int
}

func (bf BoxFeatures) countFeatures() (boxFeatureCount, error) {
	bfc := boxFeatureCount{}
	if len(bf.sx) != len(bf.sX) {
		return bfc, fmt.Errorf("X surface points not matched. %d minor, %d major", len(bf.sx), len(bf.sX))
	}
	if len(bf.sy) != len(bf.sY) {
		return bfc, fmt.Errorf("Ysurface points not matched. %d minor, %d major", len(bf.sy), len(bf.sY))
	}
	if len(bf.sz) != len(bf.sZ) {
		return bfc, fmt.Errorf("Z surface points not matched. %d minor, %d major", len(bf.sz), len(bf.sZ))
	}
	if len(bf.exZ) != len(bf.eXz) || len(bf.exz) != len(bf.eXZ) {
		return bfc, errors.New("XZ edge points not matched")
	}
	if len(bf.exy) != len(bf.eXY) || len(bf.exY) != len(bf.eXy) {
		return bfc, errors.New("XY edge points not matched")
	}
	if len(bf.eyz) != len(bf.eYZ) || len(bf.eyZ) != len(bf.eYz) {
		return bfc, errors.New("YZ edge points not matched")
	}
	bfc.sx = len(bf.sx) - len(bf.eyz) - len(bf.eYz) - len(bf.eyZ) - len(bf.eYZ) + 4
	bfc.sy = len(bf.sy) - len(bf.exz) - len(bf.eXz) - len(bf.exZ) - len(bf.eXZ) + 4
	bfc.sz = len(bf.sz) - len(bf.exy) - len(bf.eXy) - len(bf.exY) - len(bf.eXY) + 4
	bfc.ne1 = len(bf.exZ) - 2
	bfc.ne2 = len(bf.exz) - 2
	bfc.ne3 = len(bf.exy) - 2
	bfc.ne4 = len(bf.exY) - 2
	bfc.ne5 = len(bf.eyZ) - 2
	bfc.ne6 = len(bf.eyz) - 2
	return bfc, nil
}

func (bf BoxFeatures) expectedEquations() (int, error) {
	bfc, err := bf.countFeatures()
	if err != nil {
		return 0, err
	}
	return bfc.sx + bfc.sy + bfc.sz + bfc.ne1 + bfc.ne2 + bfc.ne3 + bfc.ne4 + bfc.ne5 + bfc.ne6 + 4, nil
}

// PeriodicBoundaryConditions generates Lagrange Peridic boundary conditions.
// Returned Cn SparseAccum represents a nEquations by NumberOfNodes*dofs.Count() matrix.
func (bf BoxFeatures) PeriodicBoundaryConditions(dofs fem.DofsFlag) (Cn expmat.SparseAccum, nEquations int, err error) {
	if dofs != fem.DofPos {
		return Cn, 0, errors.New("currently only support full displacement dofs for periodic boundary conditions")
	}
	rowsExpected, err := bf.expectedEquations()
	if err != nil {
		return Cn, 0, err
	}

	// Constraint matrix will be of size nEquations * NumberOfDofsInModel
	sz := rowsExpected * len(bf.nodes) * dofs.Count()
	Cn = expmat.NewSparseAccum(sz)

	// Constrain surfaces.
	var nsx, nsy, nsz int
	currentRow := 0
	for _, ix := range bf.sx {
		p := bf.nodes[ix]
		for _, iX := range bf.sX {
			P := bf.nodes[iX]
			if p.Z == P.Z && p.Y == P.Y && bf.bounds.Min.Z < p.Z && p.Z < bf.bounds.Max.Z &&
				bf.bounds.Min.Y < p.Y && p.Y < bf.bounds.Max.Y {
				linkPeriodicDisplacements(Cn, currentRow, ix, iX)
				currentRow++
				nsx++
				break
			}
			if iX == len(bf.sX)-1 {
				return Cn, 0, fmt.Errorf("minor X surface node %d has no major pair", ix)
			}
		}
	}
	// Y Surface displacement constraint.
	for _, iy := range bf.sy {
		p := bf.nodes[iy]
		for _, iY := range bf.sY {
			P := bf.nodes[iY]
			if p.Z == P.Z && p.X == P.X && bf.bounds.Min.Z < p.Z && p.Z < bf.bounds.Max.Z &&
				bf.bounds.Min.X < p.X && p.X < bf.bounds.Max.X {
				linkPeriodicDisplacements(Cn, currentRow, iy, iY)
				currentRow++
				nsy++
				break
			}
			if iY == len(bf.sY)-1 {
				return Cn, 0, fmt.Errorf("minor Y surface node %d has no major pair", iy)
			}
		}
	}
	// Z Surface displacement constraint.
	for _, iz := range bf.sz {
		p := bf.nodes[iz]
		for _, iZ := range bf.sZ {
			P := bf.nodes[iZ]
			if p.Y == P.Y && p.X == P.X && bf.bounds.Min.X < p.X && p.X < bf.bounds.Max.X &&
				bf.bounds.Min.Y < p.Y && p.Y < bf.bounds.Max.Y {
				linkPeriodicDisplacements(Cn, currentRow, iz, iZ)
				currentRow++
				nsz++
				break
			}
			if iZ == len(bf.sZ)-1 {
				return Cn, 0, fmt.Errorf("minor Z surface node %d has no major pair", iz)
			}
		}
	}

	// Constrain edges.
	var dim = func(r rune) int { return int(r - 'X') } // returns 0,1,2 with arguments 'X', 'Y' and 'Z'
	ne1 := linkRUCEdge(Cn, bf.nodes, bf.exZ, bf.eXz, currentRow, dim('Y'), bf.bounds)
	currentRow += ne1
	ne2 := linkRUCEdge(Cn, bf.nodes, bf.exz, bf.eXZ, currentRow, dim('Y'), bf.bounds)
	currentRow += ne2
	ne3 := linkRUCEdge(Cn, bf.nodes, bf.exy, bf.eXY, currentRow, dim('Z'), bf.bounds)
	currentRow += ne3
	ne4 := linkRUCEdge(Cn, bf.nodes, bf.exY, bf.eXy, currentRow, dim('Z'), bf.bounds)
	currentRow += ne4
	ne5 := linkRUCEdge(Cn, bf.nodes, bf.eyZ, bf.eYz, currentRow, dim('X'), bf.bounds)
	currentRow += ne5
	ne6 := linkRUCEdge(Cn, bf.nodes, bf.eyz, bf.eYZ, currentRow, dim('X'), bf.bounds)
	currentRow += ne6
	// Constrain corners.
	linkPeriodicDisplacements(Cn, currentRow, bf.cxyZ, bf.cXYz)
	currentRow++
	linkPeriodicDisplacements(Cn, currentRow, bf.cxyz, bf.cXYZ)
	currentRow++
	linkPeriodicDisplacements(Cn, currentRow, bf.cXyz, bf.cxYZ)
	currentRow++
	linkPeriodicDisplacements(Cn, currentRow, bf.cxYz, bf.cXyZ)
	currentRow++
	if currentRow != rowsExpected {
		return Cn, currentRow * dofs.Count(), fmt.Errorf("expected %d equations, got %d", rowsExpected, currentRow)
	}
	return Cn, rowsExpected * dofs.Count(), nil
}

func linkRUCEdge(NN expmat.SparseAccum, nodes []r3.Vec, e1, e2 []int, rows, crossDim int, modelSize r3.Box) (proc int) {
	if crossDim < 0 || crossDim > 2 {
		panic("bad cross dimension where 0,1,2 corresponds to X,Y,Z")
	}
	const (
		VecSize   = unsafe.Sizeof(r3.Vec{})
		VecOffset = unsafe.Alignof(r3.Vec{}.X)
	)
	rowsStart := rows
	nodePtr := uintptr(unsafe.Pointer(&nodes[0]))
	offset := uintptr(VecOffset * uintptr(crossDim))
	maxPtr := uintptr(unsafe.Pointer(&modelSize.Max))
	maxDim := *(*float64)(unsafe.Pointer(maxPtr + offset))
	minPtr := uintptr(unsafe.Pointer(&modelSize.Min))
	minDim := *(*float64)(unsafe.Pointer(minPtr + offset))
	for _, i1 := range e1 {
		// very unsafe. very sharp.
		pdim := *(*float64)(unsafe.Pointer(nodePtr + VecSize*uintptr(i1) + offset))
		for _, i2 := range e2 {
			Pdim := *(*float64)(unsafe.Pointer(nodePtr + VecSize*uintptr(i2) + offset))
			if pdim == Pdim && minDim < pdim && pdim < maxDim {
				linkPeriodicDisplacements(NN, rows, i1, i2)
				rows++
			}
		}
	}
	return rows - rowsStart
}

func linkPeriodicDisplacements(NN expmat.SparseAccum, r, i1, i2 int) {
	r3 := r * 3
	i13 := i1 * 3
	i23 := i2 * 3
	offset := r * 6
	NN.Set(offset, r3, i13, -1)
	NN.Set(offset+1, r3, i23, 1)
	NN.Set(offset+2, r3+1, i13+1, -1)
	NN.Set(offset+3, r3+1, i23+1, 1)
	NN.Set(offset+4, r3+2, i13+2, -1)
	NN.Set(offset+5, r3+2, i23+2, 1)
}

func NewBoxFeatures(nodes []r3.Vec, modelSize r3.Box) (bf BoxFeatures) {
	bf.nodes = nodes
	bf.bounds = modelSize
	for i, n := range nodes {
		// Sure this loop is ugly, but it should consume less energy
		// than having multiple loops since less compares. Save the trees?
		xeq0 := n.X == modelSize.Min.X
		yeq0 := n.Y == modelSize.Min.Y
		zeq0 := n.Z == modelSize.Min.Z
		xeqL := n.X == modelSize.Max.X
		yeqL := n.Y == modelSize.Max.Y
		zeqL := n.Z == modelSize.Max.Z
		if xeq0 {
			bf.sx = append(bf.sx, i)
			if yeq0 {
				bf.exy = append(bf.exy, i)
			} else if yeqL {
				bf.exY = append(bf.exY, i)
			}
			if zeq0 {
				bf.exz = append(bf.exz, i)
			} else if zeqL {
				bf.exZ = append(bf.exZ, i)
			}
		}
		if xeqL {
			bf.sX = append(bf.sX, i)
		}
		if yeq0 {
			bf.sy = append(bf.sy, i)
			if xeqL {
				bf.eXy = append(bf.eXy, i)
			}
			if zeq0 {
				bf.eyz = append(bf.eyz, i)
			} else if zeqL {
				bf.eyZ = append(bf.eyZ, i)
			}
		}
		if yeqL {
			bf.sY = append(bf.sY, i)
			if xeqL {
				bf.eXY = append(bf.eXY, i)
			}
			if zeq0 {
				bf.eYz = append(bf.eYz, i)
			} else if zeqL {
				bf.eYZ = append(bf.eYZ, i)
			}
		}
		if zeq0 {
			bf.sz = append(bf.sz, i)
			if xeqL {
				bf.eXz = append(bf.eXz, i)
			}
			if yeq0 && xeq0 {
				bf.cxyz = i
			} else if yeq0 && xeqL {
				bf.cXyz = i
			} else if yeqL && xeq0 {
				bf.cxYz = i
			} else if yeqL && xeqL {
				bf.cXYz = i
			}
		}
		if zeqL {
			bf.sZ = append(bf.sZ, i)
			if xeqL {
				bf.eXZ = append(bf.eXZ, i)
			}
			if yeq0 && xeq0 {
				bf.cxyZ = i
			} else if yeq0 && xeqL {
				bf.cXyZ = i
			} else if yeqL && xeq0 {
				bf.cxYZ = i
			} else if yeqL && xeqL {
				bf.cXYZ = i
			}
		}
	}
	return bf
}
