package solids_test

import (
	"math/rand"
	"testing"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/constitution/solids"
	"github.com/soypat/go-fem/elements"
	"github.com/soypat/lap"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

func TestPatchAxisymmetric(t *testing.T) {
	const (
		// Washer dimensions. ri=internal radius; lx=washer radial width; ly=washer height.
		ri, lx, ly = 1, 0.24, 0.12
	)
	rng := rand.New(rand.NewSource(1))
	randf := func() float64 { return rng.Float64() }

	nodes := []r3.Vec{
		{X: ri, Y: 0},
		{X: ri + lx, Y: 0},
		{X: ri + lx, Y: ly},
		{X: ri, Y: ly},
		// Internal nodes for the washer.
		{X: ri + lx/3 + (randf()-0.5)*lx/4, Y: ly/3 + (randf()-0.5)*ly/4},
		{X: ri + lx*2/3 - (randf()-0.5)*lx/4, Y: ly/3 + (randf()-0.5)*ly/4},
		{X: ri + lx*2/3 - (randf()-0.5)*lx/4, Y: ly*2/3 - (randf()-0.5)*ly/4},
		{X: ri + lx/3 + (randf()-0.5)*lx/4, Y: ly*2/3 - (randf()-0.5)*ly/4},
	}

	// Create a patch with 4 nodes and 2 elements.
	q4elems := [][4]int{
		{0, 1, 5, 4},
		{1, 2, 6, 5},
		{6, 2, 3, 7},
		{0, 4, 7, 3},
		{4, 5, 6, 7},
	}
	// getelement := func(i int) (elem []int) { return q4elems[i][:] }
	getelementR3 := func(i int) (elem []int, xC r3.Vec, yC r3.Vec) { return q4elems[i][:], r3.Vec{}, r3.Vec{} }
	material := solids.Isotropic{E: 1e6, Poisson: 0.25}
	elemtype := elements.Quad4{}
	ga := fem.NewGeneralAssembler(nodes, elemtype.Dofs())
	err := ga.AddIsoparametric(elemtype, material.Axisymmetric(), len(q4elems), getelementR3)
	if err != nil {
		t.Fatal(err)
	}
	displacements := lap.NewDenseVector(ga.TotalDofs(), nil)
	fix := fem.NewFixity(fem.DofPosX|fem.DofPosY, len(nodes))
	// Displacement loads for case 1.
	for i, node := range nodes {
		isTopNode := node.Y == ly
		isBottomNode := node.Y == 0
		isLeftNode := node.X == ri
		isRightNode := node.X == ri+lx
		isExternalNode := isTopNode || isBottomNode || isLeftNode || isRightNode
		if isExternalNode {
			displacements.SetVec(i*2, node.X/1000)
			displacements.SetVec(i*2+1, (node.Y+node.X)/1000)
			fix.Fix(i, fem.DofPosX|fem.DofPosY)
		}
	}
	free := fix.FreeDofs()
	fixed := fix.FixedDofs()
	loadsFromDispl := lap.NewDenseVector(len(free), nil)
	loadsFromDispl.MulVec(lap.Slice(ga.Ksolid(), free, fixed), lap.SliceVec(displacements, fixed))

	// 0 imposed initial imposed loads.  K(free,free)\(R(free)-K(free,fix)*D(fix))
	imposedLoads := lap.NewDenseVector(len(free), nil)
	imposedLoads.SubVec(imposedLoads, loadsFromDispl)

	resultDisplacements := mat.NewVecDense(len(free), nil)
	resultDisplacements.SolveVec(lapmat{lap.Slice(ga.Ksolid(), free, free)}, lapvec{imposedLoads})
	// Extract the solution to obtain all displacements, both fixed and free.
	lap.SliceVec(displacements, free).CopyVec(resultDisplacements)
	// Compute Reaction forces.
	forces := lap.NewDenseVector(ga.TotalDofs(), nil)
	forces.MulVec(ga.Ksolid(), displacements)
}

type lapvec struct {
	lap.Vector
}

func (v lapvec) T() mat.Matrix {
	return lapmat{lap.T(v.Vector)}
}

type lapmat struct {
	lap.Matrix
}

func (m lapmat) T() mat.Matrix {
	if ter, ok := m.Matrix.(mat.Matrix); ok {
		return ter.T()
	}
	return lapmat{lap.T(m.Matrix)}
}
