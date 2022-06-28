package tools

import (
	"errors"
	"open-anno/pkg/io"

	"github.com/go-playground/validator/v10"
)

type Any2AnyParam struct {
	Input string `validate:"required"`
	Ouput string `validate:"required"`
}

func (this Any2AnyParam) Valid() error {
	validate := validator.New()
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func Run() error {
	return errors.New("Not Impl")
}

type Av2AnyParam struct {
	Any2AnyParam
}

func (this Av2AnyParam) Variants() (io.Variants, error) {
	variants := make(io.Variants, 0)
	reader, err := io.NewIoReader(this.Input)
	if err != nil {
		return variants, err
	}
	defer reader.Close()
	scanner := io.NewVarScanner(reader)
	for scanner.Scan() {
		row, err := scanner.Row()
		if err != nil {
			return variants, err
		}
		variants = append(variants, row)
	}
	return variants, err
}
