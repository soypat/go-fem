package main

import (
	_ "embed"
	"errors"
	"fmt"
	"io"
	"strings"

	"gonum.org/v1/gonum/spatial/r3"
)

func main() {
	nodes, q8 := feaModel()
	fmt.Println("nodes:", len(nodes), "elements:", len(q8))
}

var (
	//go:embed bimetallic_nod.tsv
	_bimetallicNodes string
	//go:embed bimetallic_elem.tsv
	_bimetallicElem string
)

func feaModel() (nodes []r3.Vec, q8 [][8]int) {
	// preallocate reasonable size for warm start to appends.
	nodes = make([]r3.Vec, 0, 512)
	q8 = make([][8]int, 0, 256)
	r := strings.NewReader(_bimetallicNodes)
	// nodemap maps the model's node number to the index in the nodes slice.
	var nodemap = make(map[int]int)
	var err error
	for err == nil {
		var x, y, z float64
		var n, ignore int
		_, err = fmt.Fscanf(r, "%d\t%f\t%f\t%f\t%d\n", &n, &x, &y, &z, &ignore)
		if err == nil {
			nodemap[n] = len(nodes)
			nodes = append(nodes, r3.Vec{X: x, Y: y, Z: z})
		}
	}
	err = nil
	r = strings.NewReader(_bimetallicElem)
	for err == nil {
		var n, n1, n2, n3, n4, n5, n6, n7, n8, ignore int
		_, err = fmt.Fscanf(r, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", &n, &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &ignore)
		if err == nil {
			q8 = append(q8, [8]int{nodemap[n1], nodemap[n2], nodemap[n3], nodemap[n4], nodemap[n5], nodemap[n6], nodemap[n7], nodemap[n8]})
		} else if !errors.Is(err, io.EOF) {
			panic(err)
		}
	}
	return nodes, q8
}
