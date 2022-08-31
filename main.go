package main

import (
	"log"
	"open-anno/cmd"
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
		Use:   "pre",
		Short: "Prepare database",
	}
	cmd.AddCommand(pre.NewPreGeneBasedCmd())
	cmd.AddCommand(pre.NewIndexDatabaseCmd())
	return cmd
}

func NewAnnoCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "anno",
		Short: "Annotate variants",
	}
	cmd.AddCommand(anno.NewAnnoCmd("snv"))
	cmd.AddCommand(anno.NewAnnoCmd("cnv"))
	cmd.AddCommand(anno.NewMergeCmd())
	return cmd
}

func NewToolsCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "tools",
		Short: "Tools",
	}
	cmd.AddCommand(tools.NewVCF2AnnoInputCmd())
	cmd.AddCommand(tools.NewAnnoInput2VCFCmd())
	return cmd
}

func init() {
	RootCmd = NewRootCmd()
	RootCmd.AddCommand(NewPreCmd())
	RootCmd.AddCommand(NewAnnoCmd())
	RootCmd.AddCommand(NewToolsCmd())
	RootCmd.AddCommand(cmd.NewTestCmd())
}

func main() {
	err := RootCmd.Execute()
	if err != nil {
		log.Panic(err)
	}
}
