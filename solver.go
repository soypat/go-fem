package fem

import (
	"errors"
	"fmt"

	"github.com/soypat/go-fem/exp/expmat"
	"gonum.org/v1/gonum/mat"
)

type Solver interface {
	mat.Vector
	SolveVec(mat.Matrix, mat.Vector) error
}

type GeneralSolution struct {
	ebc      EssentialBC
	solution mat.Vector
}

// Solve solves a finite element problem given global stiffness matrix and
// global boundary conditions. Description of arguments:
//  - A linear system solver s
//  - The global stiffness matrix Kglobal
//  - Natural or Neumann boundary conditions as globalNaturalBC. These are loads
//  in solid mechanics, heat flow in thermal problems or current in electrostatic problems.
//  - Essential or Dirichlet boundary conditions as eBC. These correspond to
//  imposed displacements in solid mechanics or fixed temperatures in thermal problems.
func Solve(s Solver, Kglobal mat.Matrix, globalNaturalBC mat.Vector, eBC EssentialBC) (*GeneralSolution, error) {
	if s == nil || Kglobal == nil || globalNaturalBC == nil || eBC == nil {
		panic("nil argument to Solve")
	}
	r, c := Kglobal.Dims()
	if r == 0 || c == 0 {
		return nil, errors.New("got zero dimension in Kglobal")
	}
	if r != c {
		return nil, fmt.Errorf("Kglobal not square. got %d by $d", r, c)
	}
	if globalNaturalBC.Len() != r {
		return nil, fmt.Errorf("Kglobal dimension %d should be equal to globalNaturalBC length %d", r, globalNaturalBC.Len())
	}
	dofs := eBC.Dofs()
	dofsPerNode := dofs.Count()
	totalDofs := dofsPerNode * eBC.NumberOfNodes()
	if totalDofs != r {
		return nil, fmt.Errorf("Kglobal dimension %d should be equal to amount of dofs in EssentialBC %d", r, totalDofs)
	}
	// fixed will contain essentialBC information since we do not support imposed essentialBCs.
	fixed := make([]bool, totalDofs)
	err := setBoolEssentialBC(fixed, eBC)
	if err != nil {
		return nil, err
	}
	err = s.SolveVec(
		expmat.BoolIndexed(Kglobal, true, fixed, fixed),
		expmat.BoolIndexedVec(globalNaturalBC, true, fixed),
	)
	if err != nil {
		return nil, fmt.Errorf("Solver returned an error during Solve: %w", err)
	}
	return &GeneralSolution{
		solution: s,
		ebc:      eBC,
	}, nil
}

func (gs *GeneralSolution) DenseSolution() *mat.VecDense {
	dofsPerNode := gs.ebc.Dofs().Count()
	totalDofs := dofsPerNode * gs.ebc.NumberOfNodes()
	fixed := make([]bool, totalDofs)
	err := setBoolEssentialBC(fixed, gs.ebc)
	if err != nil {
		panic(err)
	}
	solution := mat.NewVecDense(totalDofs, nil)
	expmat.BoolSetVec(solution, gs.solution, true, fixed)
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
