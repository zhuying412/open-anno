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
	cmd.AddCommand(pre.NewPreGenePredCmd())
	cmd.AddCommand(pre.NewIdxDBCmd())
	cmd.AddCommand(pre.NewGeneCmd())
	cmd.AddCommand(pre.NewPreClinvarCmd())
	cmd.AddCommand(pre.NewPre1000GenomesCmd())
	return cmd
}

func NewAnnoCmd() *cobra.Command {
	annoGB := &cobra.Command{
		Use:   "gb",
		Short: "Annotate GeneBased",
	}
	annoGB.AddCommand(anno.NewAnnoGBSnvCmd())
	annoGB.AddCommand(anno.NewAnnoGBCnvCmd())
	annoBatch := &cobra.Command{
		Use:   "batch",
		Short: "Annotate Batch",
	}
	annoBatch.AddCommand(anno.NewAnnoBatchSnvCmd())
	annoBatch.AddCommand(anno.NewAnnoBatchCnvCmd())
	cmd := &cobra.Command{
		Use:   "anno",
		Short: "Annotate variants",
	}
	cmd.AddCommand(annoGB)
	cmd.AddCommand(anno.NewAnnoFBCmd())
	cmd.AddCommand(anno.NewAnnoRBCmd())
	cmd.AddCommand(annoBatch)
	return cmd
}

func NewToolsCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "tools",
		Short: "Tools",
	}
	cmd.AddCommand(tools.NewAggsCmd())
	cmd.AddCommand(tools.NewMergeCmd())
	cmd.AddCommand(tools.NewRepTransCmd())
	cmd.AddCommand(tools.NewExonBedCmd())
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
