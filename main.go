package main

import (
	"log"
	"open-anno/cmd/anno"
	"open-anno/cmd/pre"

	"github.com/spf13/cobra"
)

var RootCmd *cobra.Command

func NewRootCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "OpenAnno",
		Short: "A Genome Annotate Tool",
	}
	return cmd
}

func NewPreCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "Pre",
		Short: "Prepare database",
	}
	cmd.AddCommand(pre.NewPreGeneBasedCmd())
	return cmd
}

func NewAnnoCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "Anno",
		Short: "Annotate variants",
	}
	cmd.AddCommand(anno.NewAnnoGeneBasedCmd("snv"))
	cmd.AddCommand(anno.NewAnnoGeneBasedCmd("cnv"))
	return cmd
}

func init() {
	RootCmd = NewRootCmd()
	RootCmd.AddCommand(NewPreCmd())
	RootCmd.AddCommand(NewAnnoCmd())
}

func main() {
	err := RootCmd.Execute()
	if err != nil {
		log.Panic(err)
	}
}
