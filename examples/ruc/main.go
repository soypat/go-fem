package main

import (
	_ "embed"
	"encoding/csv"
	"fmt"
	"io"
	"strconv"
	"strings"
	"time"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/elements"
	"gonum.org/v1/gonum/spatial/r3"
)

func main() {
	nodes, elem := feaModel()
	ga := fem.NewGeneralAssembler(nodes, fem.DofU)
	tstart := time.Now()
	err := ga.AddIsoparametric3(elements.Hexa8{}, fem.IsotropicMaterial{E: 4.8e3, Poisson: 0.34}, len(elem), func(i int) ([]int, r3.Vec, r3.Vec) {
		return elem[i][:], r3.Vec{}, r3.Vec{}
	})

	if err != nil {
		panic(err)
	}
	dofs, _ := ga.Ksolid().Dims()
	fmt.Printf("assembled %d dofs in %s", dofs, time.Since(tstart))
}

var (
	//go:embed crettonRUC_70_elements.tsv
	_rucElem string
	//go:embed crettonRUC_70_nodes.tsv
	_rucNodes string
)

func feaModel() (nodes []r3.Vec, h8 [][8]int) {
	// preallocate reasonable size for warm start to appends.
	nodes = make([]r3.Vec, 0, 512)
	h8 = make([][8]int, 0, 256)
	r := csv.NewReader(strings.NewReader(_rucNodes))
	r.Comma = '\t'
	r.ReuseRecord = true
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}

		x, err := strconv.ParseFloat(record[1], 64)
		if err != nil {
			panic(err)
		}
		y, _ := strconv.ParseFloat(record[2], 64)
		z, _ := strconv.ParseFloat(record[3], 64)
		nodes = append(nodes, r3.Vec{X: x, Y: y, Z: z})
	}
	r = csv.NewReader(strings.NewReader(_rucElem))
	r.Comma = '\t'
	r.ReuseRecord = true
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		var elem [8]int
		for i := 0; i < 8; i++ {
			elem[i], err = strconv.Atoi(record[i+1])
			if err != nil {
				panic(err)
			}
			elem[i]-- // to account for 1 indexing.
		}
		h8 = append(h8, elem)
	}
	return nodes, h8
}
