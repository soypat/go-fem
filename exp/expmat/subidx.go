package expmat

import "gonum.org/v1/gonum/mat"

// SubIdx allows
type SubIdx struct {
	cidx, ridx []int
	m          mat.Matrix
}

func (bm *SubIdx) At(i, j int) float64 { return bm.m.At(bm.ridx[i], bm.cidx[j]) }
func (bm *SubIdx) Dims() (int, int)    { return len(bm.ridx), len(bm.cidx) }
func (bm *SubIdx) T() mat.Matrix {
	return mat.Transpose{Matrix: bm}
}

type subIdxVec struct {
	SubIdx
}

func (bm *subIdxVec) Len() int            { return len(bm.ridx) }
func (bm *subIdxVec) AtVec(i int) float64 { return bm.At(i, 0) }

func NewSubIdxVec(ridx []int, m mat.Matrix) mat.Vector {
	r, c := m.Dims()
	if r == 1 {
		return &subIdxVec{*NewSubIdx([]int{0}, ridx, m)}
	} else if c != 1 {
		panic("matrix should be 1xn or nx1 to create a Vector")
	}
	return &subIdxVec{*NewSubIdx(ridx, []int{0}, m)}
}

func NewSubIdx(ridx, cidx []int, m mat.Matrix) *SubIdx {
	if len(ridx) == 0 || len(cidx) == 0 {
		panic("cannot have zero dimension SubIdx")
	}
	r, c := m.Dims()
	for _, ir := range ridx {
		if ir >= r {
			panic(mat.ErrRowAccess)
		}
	}
	for _, ic := range cidx {
		if ic >= c {
			panic(mat.ErrColAccess)
		}
	}
	return &SubIdx{ridx: ridx, cidx: cidx, m: m}
}

func BoolIndexedVec(m mat.Vector, inv bool, br []bool) mat.Vector {
	r, c := m.Dims()
	if c != 1 && r != 1 {
		panic("not a vector")
	}
	if c != 1 {
		return &subIdxVec{*BoolIndexed(m.T(), inv, br, []bool{!inv})}
	}
	return &subIdxVec{*BoolIndexed(m, inv, br, []bool{!inv})}
}

func BoolSetVec(dst *mat.VecDense, src mat.Vector, inv bool, br []bool) {
	if len(br) != dst.Len() {
		panic("bad []bool len. must match dst")
	}
	sum := 0
	for i := range br {
		if br[i] != inv {
			sum++
		}
	}
	if sum != src.Len() {
		panic("amount of true values in br must match length of src")
	}
	sum = 0
	for i := range br {
		if br[i] != inv {
			dst.SetVec(i, src.AtVec(sum))
			sum++
		}
	}
}

func BoolIndexed(m mat.Matrix, inv bool, br, bc []bool) *SubIdx {
	r, c := m.Dims()
	if len(br) != r || len(bc) != c {
		panic("bad dim")
	}
	sm := SubIdx{
		ridx: make([]int, 0, r/16),
		cidx: make([]int, 0, c/16),
		m:    m,
	}
	for i, b := range br {
		if b != inv {
			sm.ridx = append(sm.ridx, i)
		}
	}
	for i, b := range bc {
		if b != inv {
			sm.cidx = append(sm.cidx, i)
		}
	}
	return &sm
}
