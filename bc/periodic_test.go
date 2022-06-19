package bc_test

import (
	"testing"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/bc"
	"github.com/soypat/go-fem/elements"
	"github.com/soypat/go-fem/exp/expmat"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

func TestAssemblerPeriodicBoundaryConditions(t *testing.T) {
	var (
		modelBounds = r3.Box{Max: r3.Vec{X: 10, Y: 10, Z: 10}}
		nodes       = rucNodes
		material    = fem.IsotropicMaterial{E: 14e3, Poisson: 0.2}
		elemType    = elements.Hexa8{}
	)
	ga := fem.NewGeneralAssembler(nodes, fem.DofPos)
	err := ga.AddIsoparametric3(elemType, material, len(rucElements), func(i int) (elem []int, xC r3.Vec, yC r3.Vec) {
		return rucElements[i][:], xC, yC
	})
	if err != nil {
		t.Fatal(err)
	}
	bf := bc.NewBoxFeatures(nodes, modelBounds)
	Cn, nEquations, err := bf.PeriodicBoundaryConditions(fem.DofPos)
	if err != nil {
		t.Fatal(err)
	}
	Kglobal := ga.Ksolid()
	r, c := Kglobal.Dims()
	Kglobal.Resize(r+nEquations, c+nEquations)
	Kglobal.GeneralAccumulate(false, r, 0, Cn)
	Kglobal.GeneralAccumulate(true, 0, c, Cn)
	ndofs, _ := Kglobal.Dims()
	fixedDofs := make([]bool, ndofs)
	boolsAreFixed := true
	// Fix one node of the model so it cannot translate as a rigid body.
	fixedDofs[len(nodes)*3-1] = true
	fixedDofs[len(nodes)*3-2] = true
	fixedDofs[len(nodes)*3-3] = true
	loads := make([]float64, ndofs)
	imposedLoads := loads[len(loads)-nEquations:]
	// Cruc := mat.NewDense(6, 6, nil)
	displacements := mat.NewVecDense(ndofs, nil)
	for rucCase := 0; rucCase < 6; rucCase++ {
		err = bf.PeriodicImposedDisplacements(imposedLoads, rucCase, 0.1)
		if err != nil {
			t.Fatal(err)
		}
		// Solve system.
		var freeDisplacements mat.VecDense
		loadVec := mat.NewVecDense(ndofs, loads)
		freeDisplacements.SolveVec(
			expmat.BoolIndexed(Kglobal, boolsAreFixed, fixedDofs, fixedDofs),
			expmat.BoolIndexedVec(loadVec, boolsAreFixed, fixedDofs),
		)
		expmat.BoolSetVec(displacements, &freeDisplacements, boolsAreFixed, fixedDofs)

		// Extract stresses.
		// unod := elemType.IsoparametricNodes()
	}

}

var rucElements = [][8]int{
	0:  {0, 1, 4, 3, 9, 10, 13, 12},
	1:  {1, 2, 5, 4, 10, 11, 14, 13},
	2:  {3, 4, 7, 6, 12, 13, 16, 15},
	3:  {4, 5, 8, 7, 13, 14, 17, 16},
	4:  {9, 10, 13, 12, 18, 19, 22, 21},
	5:  {10, 11, 14, 13, 19, 20, 23, 22},
	6:  {12, 13, 16, 15, 21, 22, 25, 24},
	7:  {13, 14, 17, 16, 22, 23, 26, 25},
	8:  {27, 28, 1, 0, 30, 31, 10, 9},
	9:  {28, 29, 2, 1, 31, 32, 11, 10},
	10: {30, 31, 10, 9, 33, 34, 19, 18},
	11: {31, 32, 11, 10, 34, 35, 20, 19},
	12: {36, 37, 3, 6, 38, 39, 12, 15},
	13: {37, 27, 0, 3, 39, 30, 9, 12},
	14: {38, 39, 12, 15, 40, 41, 21, 24},
	15: {39, 30, 9, 12, 41, 33, 18, 21},
	16: {42, 43, 7, 8, 44, 45, 16, 17},
	17: {43, 36, 6, 7, 45, 38, 15, 16},
	18: {44, 45, 16, 17, 46, 47, 25, 26},
	19: {45, 38, 15, 16, 47, 40, 24, 25},
	20: {29, 48, 5, 2, 32, 49, 14, 11},
	21: {48, 42, 8, 5, 49, 44, 17, 14},
	22: {32, 49, 14, 11, 35, 50, 23, 20},
	23: {49, 44, 17, 14, 50, 46, 26, 23},
	24: {51, 52, 28, 27, 54, 55, 31, 30},
	25: {52, 53, 29, 28, 55, 56, 32, 31},
	26: {54, 55, 31, 30, 57, 58, 34, 33},
	27: {55, 56, 32, 31, 58, 59, 35, 34},
	28: {60, 61, 37, 36, 62, 63, 39, 38},
	29: {61, 51, 27, 37, 63, 54, 30, 39},
	30: {62, 63, 39, 38, 64, 65, 41, 40},
	31: {63, 54, 30, 39, 65, 57, 33, 41},
	32: {66, 67, 43, 42, 68, 69, 45, 44},
	33: {67, 60, 36, 43, 69, 62, 38, 45},
	34: {68, 69, 45, 44, 70, 71, 47, 46},
	35: {69, 62, 38, 45, 71, 64, 40, 47},
	36: {53, 72, 48, 29, 56, 73, 49, 32},
	37: {72, 66, 42, 48, 73, 68, 44, 49},
	38: {56, 73, 49, 32, 59, 74, 50, 35},
	39: {73, 68, 44, 49, 74, 70, 46, 50},
}

var rucNodes = []r3.Vec{
	0:  {X: 10, Y: 8, Z: 8},
	1:  {X: 10, Y: 5, Z: 8},
	2:  {X: 10, Y: 2, Z: 8},
	3:  {X: 10, Y: 8, Z: 5},
	4:  {X: 10, Y: 5, Z: 5},
	5:  {X: 10, Y: 2, Z: 5},
	6:  {X: 10, Y: 8, Z: 2},
	7:  {X: 10, Y: 5, Z: 2},
	8:  {X: 10, Y: 2, Z: 2},
	9:  {X: 5, Y: 8, Z: 8},
	10: {X: 5, Y: 5, Z: 8},
	11: {X: 5, Y: 2, Z: 8},
	12: {X: 5, Y: 8, Z: 5},
	13: {X: 5, Y: 5, Z: 5},
	14: {X: 5, Y: 2, Z: 5},
	15: {X: 5, Y: 8, Z: 2},
	16: {X: 5, Y: 5, Z: 2},
	17: {X: 5, Y: 2, Z: 2},
	18: {X: 0, Y: 8, Z: 8},
	19: {X: 0, Y: 5, Z: 8},
	20: {X: 0, Y: 2, Z: 8},
	21: {X: 0, Y: 8, Z: 5},
	22: {X: 0, Y: 5, Z: 5},
	23: {X: 0, Y: 2, Z: 5},
	24: {X: 0, Y: 8, Z: 2},
	25: {X: 0, Y: 5, Z: 2},
	26: {X: 0, Y: 2, Z: 2},
	27: {X: 10, Y: 8.337790589, Z: 8.337790589},
	28: {X: 10, Y: 5, Z: 9.72034871932508},
	29: {X: 10, Y: 1.662209411, Z: 8.337790589},
	30: {X: 5, Y: 8.337790589, Z: 8.337790589},
	31: {X: 5, Y: 5, Z: 9.72034871932508},
	32: {X: 5, Y: 1.662209411, Z: 8.337790589},
	33: {X: 0, Y: 8.337790589, Z: 8.337790589},
	34: {X: 0, Y: 5, Z: 9.72034871932508},
	35: {X: 0, Y: 1.662209411, Z: 8.337790589},
	36: {X: 10, Y: 8.337790589, Z: 1.662209411},
	37: {X: 10, Y: 9.72034871932508, Z: 5},
	38: {X: 5, Y: 8.337790589, Z: 1.662209411},
	39: {X: 5, Y: 9.72034871932508, Z: 5},
	40: {X: 0, Y: 8.337790589, Z: 1.662209411},
	41: {X: 0, Y: 9.72034871932508, Z: 5},
	42: {X: 10, Y: 1.662209411, Z: 1.662209411},
	43: {X: 10, Y: 5, Z: 0.279651280674918},
	44: {X: 5, Y: 1.662209411, Z: 1.662209411},
	45: {X: 5, Y: 5, Z: 0.279651280674918},
	46: {X: 0, Y: 1.662209411, Z: 1.662209411},
	47: {X: 0, Y: 5, Z: 0.279651280674918},
	48: {X: 10, Y: 0.279651280674919, Z: 5},
	49: {X: 5, Y: 0.279651280674919, Z: 5},
	50: {X: 0, Y: 0.279651280674919, Z: 5},
	51: {X: 10, Y: 10, Z: 10},
	52: {X: 10, Y: 5, Z: 10},
	53: {X: 10, Y: 0, Z: 10},
	54: {X: 5, Y: 10, Z: 10},
	55: {X: 5, Y: 5, Z: 10},
	56: {X: 5, Y: 0, Z: 10},
	57: {X: 0, Y: 10, Z: 10},
	58: {X: 0, Y: 5, Z: 10},
	59: {X: 0, Y: 0, Z: 10},
	60: {X: 10, Y: 10, Z: 0},
	61: {X: 10, Y: 10, Z: 5},
	62: {X: 5, Y: 10, Z: 0},
	63: {X: 5, Y: 10, Z: 5},
	64: {X: 0, Y: 10, Z: 0},
	65: {X: 0, Y: 10, Z: 5},
	66: {X: 10, Y: 0, Z: 0},
	67: {X: 10, Y: 5, Z: 0},
	68: {X: 5, Y: 0, Z: 0},
	69: {X: 5, Y: 5, Z: 0},
	70: {X: 0, Y: 0, Z: 0},
	71: {X: 0, Y: 5, Z: 0},
	72: {X: 10, Y: 0, Z: 5},
	73: {X: 5, Y: 0, Z: 5},
	74: {X: 0, Y: 0, Z: 5},
}
