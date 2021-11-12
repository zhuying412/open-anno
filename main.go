package main

import (
	"OpenAnno/command"
	"log"
)

func main() {
	err := command.RootCmd.Execute()
	if err != nil {
		log.Panic(err)
	}
}
