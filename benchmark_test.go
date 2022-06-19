package fem_test

import (
	"fmt"
	"testing"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/elements"
	"github.com/soypat/manigold/tetra"
	"gonum.org/v1/gonum/spatial/r3"
)

//Before Sparse implementation using mat.Dense:
// cpu: AMD Ryzen 5 3400G with Radeon Vega Graphics
// BenchmarkTetra4Assembly/567_dofs-8                  1494            807242 ns/op         2580667 B/op         52 allocs/op
// BenchmarkTetra4Assembly/3723_dofs-8                   56          19811755 ns/op        110896956 B/op        56 allocs/op
// BenchmarkTetra4Assembly/27027_dofs-8                   1        1300772056 ns/op        5843683984 B/op       64 allocs/op
// BenchmarkTetra4Assembly/38073_dofs-8                   1        2628710793 ns/op        11596442256 B/op              64 allocs/op
// After that you get OOMs approaching tens of GB of memory per stiffness matrix.

// cpu: AMD Ryzen 5 2400G with Radeon Vega Graphics
// BenchmarkTetra4Assembly/105_dofs,_48_elems-8         	    1188	   1010943 ns/op	  153933 B/op	     695 allocs/op
// BenchmarkTetra4Assembly/3723_dofs,_5376_elems-8      	       9	 120619633 ns/op	10730002 B/op	   69083 allocs/op
// BenchmarkTetra4Assembly/27027_dofs,_46080_elems-8    	       1	1145319400 ns/op	87443736 B/op	  588387 allocs/op
// BenchmarkTetra4Assembly/206115_dofs,_380928_elems-8  	       1	10082928000 ns/op	706745048 B/op	 4848576 allocs/op
// BenchmarkTetra4Assembly/684723_dofs,_1299456_elems-8 	       1	37995567200 ns/op	2698002184 B/op	16588366 allocs/op
func BenchmarkTetra4Assembly(b *testing.B) {
	const dim = 1.0
	box := r3.Box{Max: r3.Vec{X: dim, Y: dim, Z: dim}}
	elemT := elements.Tetra4{}
	material := fem.IsotropicMaterial{E: 200e9, Poisson: 0.3}
	for _, div := range []float64{2, 8, 16, 32, 48} {
		b.StopTimer()
		bcc := tetra.MakeBCC(box, dim/div)
		nodes, tetras := bcc.MeshTetraBCC()
		totalDofs := len(nodes) * 3
		b.StartTimer()
		b.Run(fmt.Sprintf("%d dofs, %d elems", totalDofs, len(tetras)), func(b *testing.B) {
			var ga *fem.GeneralAssembler
			for i := 0; i < b.N; i++ {
				ga = fem.NewGeneralAssembler(nodes, fem.DofPos)
				err := ga.AddIsoparametric3(elemT, material, len(tetras), func(i int) (elem []int, xC r3.Vec, yC r3.Vec) {
					return tetras[i][:], xC, yC
				})
				if err != nil {
					b.Fatal(err)
				}
			}
			// b.Logf("%d non zero entries", ga.Ksolid().NonZero())
		})
	}
}
