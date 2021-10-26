package main

import (
	"grandanno/cmd"
	"log"
)

func main() {
	err := cmd.RootCmd.Execute()
	if err != nil {
		log.Panic(err)
	}
}
