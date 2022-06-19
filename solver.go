package fem

import (
	"errors"
	"fmt"
	"math"

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
	totalDofs := eBC.NumberOfDofs()
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

// GlobalSolution returns the globally indexed solution of the finite element problem.
func (gs *GeneralSolution) GlobalSolution() mat.Vector {
	totalDofs := gs.ebc.NumberOfDofs()
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
	for i := 0; i < eBC.NumberOfDofs(); i++ {
		bc, imposed := eBC.AtDof(i)
		if imposed != 0 {
			return errors.New("Solve does not yet support imposed essential boundary conditions")
		}
		dst[i] = bc
	}
	return nil
}

type ElementSolution struct {
	NodeStrain   []float64
	NodeStress   []float64
	QuadStrain   []float64
	QuadStress   []float64
	InterpStrain []float64
	InterpStress []float64
}

func (gs *GeneralSolution) DoStrainIsoparametric3(elemT Isoparametric3, c Constituter,
	Nelem int, getElement func(i int) (elem []int, xC, yC r3.Vec), doSol func(i int, sol ElementSolution)) error {
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
	dofMapping, err := DofMapping(dofs, elemT)
	if err != nil {
		return err
	}
	upg, wpg := elemT.Quadrature()
	if len(upg) != len(wpg) {
		return fmt.Errorf("length of element quadrature positions %d and weights %d must be equal", len(upg), len(wpg))
	}
	_, dNn := evalIsoBasis(elemT, elemT.IsoparametricNodes())
	Npg, dNpg := evalIsoBasis(elemT, upg)
	solution := gs.GlobalSolution()
	esol := ElementSolution{
		NodeStrain: make([]float64, nodesPerElem*constitutiveDim),
		NodeStress: make([]float64, nodesPerElem*constitutiveDim),
		QuadStrain: make([]float64, len(upg)*constitutiveDim),
		QuadStress: make([]float64, len(upg)*constitutiveDim),
	}
	jac := r3.NewMat(nil)
	nodeStrain := mat.NewVecDense(len(esol.NodeStrain), nil)
	nodeStress := mat.NewVecDense(len(esol.NodeStress), nil)
	interpStress := mat.NewVecDense(len(esol.InterpStress), nil)
	interpStrain := mat.NewVecDense(len(esol.InterpStrain), nil)
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
		storeElemDofs(elemDofs, element, dofMapping, dofsPerElemNode)
		// Iterate over
		for inode := 0; inode < elemT.LenNodes(); inode++ {
			dNi := dNn[inode]
			jac.Mul(dNi, elemNodes)
			dJac := jac.Det()
			if math.Abs(dJac) < 1e-12 {
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
			offset := inode * constitutiveDim
			nodeStrain.SetRawVector(blas64.Vector{
				N:    constitutiveDim,
				Inc:  1,
				Data: esol.NodeStrain[offset : offset+constitutiveDim],
			})
			nodeStrain.MulVec(B, elemDisplacements)
			nodeStress.SetRawVector(blas64.Vector{
				N:    constitutiveDim,
				Inc:  1,
				Data: esol.NodeStress[offset : offset+constitutiveDim],
			})
			nodeStress.MulVec(C, nodeStrain)
		}
		nodeStresses := mat.NewDense(nodesPerElem, constitutiveDim, esol.NodeStress)
		nodeStrains := mat.NewDense(nodesPerElem, constitutiveDim, esol.NodeStrain)
		// Quadrature interpolation.
		for i := range upg {
			offset := i * constitutiveDim
			interpStrain.SetRawVector(blas64.Vector{
				N:    constitutiveDim,
				Inc:  1,
				Data: esol.InterpStrain[offset : offset+constitutiveDim],
			})
			interpStrain.MulVec(nodeStrains, Npg[i])
			interpStress.SetRawVector(blas64.Vector{
				N:    constitutiveDim,
				Inc:  1,
				Data: esol.InterpStress[offset : offset+constitutiveDim],
			})
			interpStress.MulVec(nodeStresses, Npg[i])
		}
		// Finally we calculate stress at quadrature points which is most accurate.
		for i := range upg {
			dNi := dNpg[i]
			jac.Mul(dNi, elemNodes)
			dJac := jac.Det()
			if math.Abs(dJac) < 1e-12 {
				return fmt.Errorf("zero determinant of jacobian of element #%d, Check element shape for bad aspect ratio", iele)
			}
			err := dNxyz.Solve(jac, dNi)
			if err != nil {
				return fmt.Errorf("error calculating element #%d form factor: %s", iele, err)
			}
		}
	}
	return nil
}
