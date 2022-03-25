package main

import (
	"log"
	"open-anno/cmd/anno"
	"open-anno/cmd/pre"
	"open-anno/cmd/tools"

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
	cmd.AddCommand(pre.NewPreDatabaseCmd())
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

func NewToolsCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "Tools",
		Short: "Tools",
	}
	cmd.AddCommand(tools.NewMergeCmd())
	return cmd
}

func init() {
	RootCmd = NewRootCmd()
	RootCmd.AddCommand(NewPreCmd())
	RootCmd.AddCommand(NewAnnoCmd())
	RootCmd.AddCommand(NewToolsCmd())
}

func main() {
	err := RootCmd.Execute()
	if err != nil {
		log.Panic(err)
	}
}
