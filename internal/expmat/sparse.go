package expmat

import "gonum.org/v1/gonum/mat"

// Sparse is experimental sparse matrix.
type Sparse struct {
	m    map[[2]int]float64
	r, c int
}

var _ mat.Matrix = (*Sparse)(nil)

func (s *Sparse) At(i, j int) float64 {
	if i >= s.r {
		panic(mat.ErrRowAccess)
	} else if j >= s.c {
		panic(mat.ErrColAccess)
	}
	return s.m[[2]int{i, j}]
}

func (s *Sparse) Set(i, j int, v float64) {
	if i >= s.r {
		panic(mat.ErrRowAccess)
	} else if j >= s.c {
		panic(mat.ErrColAccess)
	}
	ix := [2]int{i, j}
	if v != 0 {
		s.m[ix] = v
	} else {
		delete(s.m, ix)
	}
}

func (s *Sparse) Dims() (int, int) {
	return s.r, s.c
}

func (s *Sparse) T() mat.Matrix { return mat.Transpose{Matrix: s} }

func NewSparse(r, c int) *Sparse {
	s := Sparse{
		r: r,
		c: c,
		m: make(map[[2]int]float64),
	}
	return &s
}

func (s *Sparse) Resize(r, c int) {
	if r < s.r || c < s.c {
		for k := range s.m {
			if k[0] >= r || k[1] >= c {
				delete(s.m, k)
			}
		}
	}
	s.r = r
	s.c = c
}

func (s *Sparse) AddData(ixs, jxs []int, data []float64) {
	if len(ixs) != len(jxs) || len(data) != len(ixs) {
		panic("length of arguments must be equal")
	}
	if len(s.m) == 0 {
		s.m = make(map[[2]int]float64, len(data)/16)
	}
	for i, v := range data {
		if v == 0 {
			continue
		}
		ix := ixs[i]
		jx := jxs[i]
		if ix >= s.r {
			panic(mat.ErrRowAccess)
		} else if jx >= s.c {
			panic(mat.ErrColAccess)
		}

		idx := [2]int{ix, jx}
		s.m[idx] += v
	}
}

// ForEachNonZero iterates over all non-zero values in sparse matrix.
// No specific ordering is guaranteed.
func (s *Sparse) ForEachNonZero(f func(i, j int, v float64)) {
	for k, v := range s.m {
		f(k[0], k[1], v)
	}
}

// NonZero returns number of non-zero elements in sparse matrix.
func (s *Sparse) NonZero() int { return len(s.m) }
