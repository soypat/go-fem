package expmat

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

func CopyBlocks(dst *mat.Dense, rows, cols int, src []mat.Matrix) error {
	if len(src) != rows*cols {
		return mat.ErrShape
	}
	var tr, tc int
	for i := 0; i < rows; i++ {
		r, _ := src[i*cols].Dims()
		tr += r
	}
	for j := 0; j < cols; j++ {
		_, c := src[j].Dims()
		tc += c
	}
	dst.ReuseAs(tr, tc)

	var br int
	for i := 0; i < rows; i++ {
		var bc int
		h, _ := src[i*cols].Dims()
		for j := 0; j < cols; j++ {
			r, c := src[i*cols+j].Dims()
			if r != h {
				return fmt.Errorf("matrix at %d,%d is wrong height: %d != %d:  %w", i, j, r, h, mat.ErrShape)
			}
			if i != 0 {
				_, w := src[j].Dims()
				if c != w {
					return fmt.Errorf("matrix at %d,%d is wrong width: %d != %d:  %w", i, j, c, w, mat.ErrShape)
				}
			}
			dst.Slice(br, br+r, bc, bc+c).(*mat.Dense).Copy(src[i*cols+j])
			bc += c
		}
		br += h
	}
	return nil
}

type Eye int

func (m Eye) At(i, j int) float64 {
	if i == j {
		return 1
	}
	return 0
}
func (m Eye) Dims() (r, c int) { return int(m), int(m) }
func (m Eye) T() mat.Matrix    { return m }

type Zero struct{ r, c int }

func (m Zero) At(i, j int) float64 { return 0 }
func (m Zero) Dims() (r, c int)    { return m.r, m.c }
func (m Zero) T() mat.Matrix       { return Zero{m.c, m.r} }

func exampleDense_copyBlocks() {
	var m mat.Dense
	r := mat.NewDense(5, 3, []float64{
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
		13, 14, 15,
	})
	err := CopyBlocks(&m, 2, 3, []mat.Matrix{
		Eye(5), Zero{5, 3}, r,
		Zero{3, 5}, Eye(3), Zero{3, 3},
	})
	if err != nil {
		fmt.Println(err)
	}
	fmt.Println(mat.Formatted(&m))
}
