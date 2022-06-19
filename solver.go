package fem

import (
	"errors"
	"fmt"

	"github.com/soypat/go-fem/exp/expmat"
	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

type Solver interface {
	mat.Vector
	SolveVec(mat.Matrix, mat.Vector) error
}

// GeneralSolution stores a finite element result from a call to Solve method.
// The zero value of GeneralSolution is ready to be used with Solve.
type GeneralSolution struct {
	ebc     EssentialBC
	sol     Solver
	nodes   []r3.Vec
	strains []float64
}

// NewGeneralSolution creates a solver with nodes and a user designated solver.
func NewGeneralSolution(solver Solver, nodes []r3.Vec) *GeneralSolution {
	if solver == nil {
		panic("nil solver argument")
	}
	return &GeneralSolution{
		sol:   solver,
		nodes: nodes,
	}
}

// Solve solves a finite element problem given global stiffness matrix and
// global boundary conditions. Description of arguments:
//  - The global stiffness matrix Kglobal
//  - Natural or Neumann boundary conditions as globalNaturalBC. These are loads
//  in solid mechanics, heat flow in thermal problems or current in electrostatic problems.
//  - Essential or Dirichlet boundary conditions as eBC. These correspond to
//  imposed displacements in solid mechanics or fixed temperatures in thermal problems.
func (gs *GeneralSolution) Solve(Kglobal mat.Matrix, globalNaturalBC mat.Vector, eBC EssentialBC) error {
	if len(gs.nodes) != 0 && len(gs.nodes) != eBC.NumberOfNodes() {
		return fmt.Errorf("GeneralSolution nodes length %d must be zero or equal to number of nodes %d expected by EssentialBC", len(gs.nodes), eBC.NumberOfNodes())
	}
	if Kglobal == nil || globalNaturalBC == nil || eBC == nil {
		panic("nil argument to Solve")
	}
	r, c := Kglobal.Dims()
	if r == 0 || c == 0 {
		return errors.New("got zero dimension in Kglobal")
	}
	if r != c {
		return fmt.Errorf("Kglobal not square. got %d by %d", r, c)
	}
	if globalNaturalBC.Len() != r {
		return fmt.Errorf("Kglobal dimension %d should be equal to globalNaturalBC length %d", r, globalNaturalBC.Len())
	}
	dofs := eBC.Dofs()
	dofsPerNode := dofs.Count()
	totalDofs := dofsPerNode * eBC.NumberOfNodes()
	if totalDofs != r {
		return fmt.Errorf("Kglobal dimension %d should be equal to amount of dofs in EssentialBC %d", r, totalDofs)
	}
	// fixed will contain essentialBC information since we do not support imposed essentialBCs.
	fixed := make([]bool, totalDofs)
	err := setBoolEssentialBC(fixed, eBC)
	if err != nil {
		return err
	}
	if gs.sol == nil {
		gs.sol = &mat.VecDense{}
	}
	err = gs.sol.SolveVec(
		expmat.BoolIndexed(Kglobal, true, fixed, fixed),
		expmat.BoolIndexedVec(globalNaturalBC, true, fixed),
	)
	if err != nil {
		return fmt.Errorf("Solver returned an error during Solve: %w", err)
	}
	gs.ebc = eBC
	return nil
}

// DenseSolution returns the globally indexed solution of the finite element problem.
func (gs *GeneralSolution) DenseSolution() *mat.VecDense {
	dofsPerNode := gs.ebc.Dofs().Count()
	totalDofs := dofsPerNode * gs.ebc.NumberOfNodes()
	fixed := make([]bool, totalDofs)
	err := setBoolEssentialBC(fixed, gs.ebc)
	if err != nil {
		panic(err)
	}
	solution := mat.NewVecDense(totalDofs, nil)
	expmat.BoolSetVec(solution, gs.sol, true, fixed)
	return solution
}

func setBoolEssentialBC(dst []bool, eBC EssentialBC) error {
	dofs := eBC.Dofs()
	dofsPerNode := dofs.Count()
	for i := 0; i < eBC.NumberOfNodes(); i++ {
		nodeDofs, imposed := eBC.AtNode(i)
		if imposed != nil {
			return errors.New("Solve does not yet support imposed essential boundary conditions")
		}
		offset := dofsPerNode * i
		dofCount := 0
		for j := 0; dofCount < dofsPerNode; j++ {
			d := DofsFlag(1 << j)
			if !dofs.Has(d) {
				continue
			}
			dst[offset+dofCount] = nodeDofs.Has(d)
			dofCount++
			if j > maxDofsPerNode {
				panic("unreachable") // Remove after tests have been added.
			}
		}
	}
	return nil
}

func (gs *GeneralSolution) StoreIsoparamtric3QuadStrains(elemT Isoparametric3, c Constituter, Nelem int, getElement func(i int) (elem []int, xC, yC r3.Vec)) error {
	if len(gs.nodes) == 0 {
		return errors.New("nodes must not be zero length before calculating isoparametric strains")
	}
	C := c.Constitutive()
	constitutiveDim, _ := C.Dims()
	var (
		dofs            = gs.ebc.Dofs()
		nodesPerElem    = elemT.LenNodes()
		dofsPerElemNode = elemT.Dofs().Count()
		dofsPerElem     = nodesPerElem * dofsPerElemNode
	)
	// helper assembler.
	ga := GeneralAssembler{
		dofs:  dofs,
		nodes: gs.nodes,
	}
	dofMapping, err := ga.DofMapping(elemT)
	if err != nil {
		return err
	}
	upg, wpg := elemT.Quadrature()
	if len(upg) != len(wpg) {
		return fmt.Errorf("length of element quadrature positions %d and weights %d must be equal", len(upg), len(wpg))
	}
	// Calculate form functions evaluated at integration points.
	N := make([]*mat.VecDense, len(upg))
	dN := make([]*mat.Dense, len(upg))
	for inod, nod := range upg {
		N[inod] = mat.NewVecDense(1, elemT.Basis(nod))
		dN[inod] = mat.NewDense(dofsPerElemNode, nodesPerElem, elemT.BasisDiff(nod))
	}

	solution := gs.DenseSolution()
	var auxVec mat.VecDense
	jac := r3.NewMat(nil)
	strain := make([]float64, Nelem*len(upg)*constitutiveDim)
	elemDofs := make([]int, nodesPerElem*dofsPerElemNode)
	elemNodesBacking := make([]float64, 3*nodesPerElem)
	elemNodes := mat.NewDense(nodesPerElem, 3, elemNodesBacking)
	B := mat.NewDense(constitutiveDim, dofsPerElem, nil)
	dNxyz := mat.NewDense(3, nodesPerElem, nil)
	elemDisplacements := expmat.NewSubIdxVec(elemDofs, solution)
	for iele := 0; iele < Nelem; iele++ {
		element, xC, yC := getElement(iele)
		if len(element) != nodesPerElem {
			return fmt.Errorf("element #%d of %d nodes expected to be of %d nodes", iele, len(element), nodesPerElem)
		}
		if xC != (r3.Vec{}) || yC != (r3.Vec{}) {
			return errors.New("isoparametric element constitutive matrix arbitrary orientation not supported yet")
		}
		storeElemNode(elemNodesBacking, gs.nodes, element)
		storeElemDofs(elemDofs, element, dofMapping, dofsPerElemNode) // This modifies elemDisplacements.
		for inode := range upg {
			dNi := dN[inode]
			jac.Mul(dNi, elemNodes)
			dJac := jac.Det()
			if dJac < 0 {
				// return fmt.Errorf("negative determinant of jacobian of element #%d, Check node ordering", iele)
			} else if dJac < 1e-12 {
				return fmt.Errorf("zero determinant of jacobian of element #%d, Check element shape for bad aspect ratio", iele)
			}
			err := dNxyz.Solve(jac, dNi)
			if err != nil {
				return fmt.Errorf("error calculating element #%d form factor: %s", iele, err)
			}
			for i := 0; i < 8; i++ {
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
			auxIdx := iele*(nodesPerElem*constitutiveDim) + inode*constitutiveDim
			// strain = B*D  where D is displacements
			s := strain[auxIdx : auxIdx+6]
			auxVec.SetRawVector(blas64.Vector{N: 6, Inc: 1, Data: s})
			// To calculate stresses later on: stress = C*B*D = C*strain
			auxVec.MulVec(B, elemDisplacements)
		}
	}
	return nil
}
