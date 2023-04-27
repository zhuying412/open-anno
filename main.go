package main

import (
	"log"
	"open-anno/cmd"
	"open-anno/cmd/anno"
	"open-anno/cmd/pre"
	"open-anno/cmd/pre/clinvar"
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
	cmd.AddCommand(pre.NewGeneCmd())

	cmd.AddCommand(pre.NewPreGnomadCmd())
	cmd.AddCommand(pre.NewPreDbnsfpCmd())
	cmd.AddCommand(pre.NewSplitVCFCmd())
	cln := &cobra.Command{
		Use:   "clinvar",
		Short: "Prepare ClinVar Database",
	}
	cln.AddCommand(clinvar.NewPreClinvarCmd())
	cln.AddCommand(clinvar.NewPrePathogenicCmd())
	cln.AddCommand(clinvar.NewPreClinvarGeneCmd())
	cmd.AddCommand(cln)
	return cmd
}

func NewAnnoCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "anno",
		Short: "Annotate variants",
	}
	cmd.AddCommand(anno.NewAnnoSnvCmd())
	cmd.AddCommand(anno.NewAnnoCnvCmd())
	return cmd
}

func NewToolsCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "tools",
		Short: "Tools",
	}
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
	// var rlim syscall.Rlimit
	// err := syscall.Getrlimit(syscall.RLIMIT_NOFILE, &rlim)
	// if err != nil {
	// 	log.Fatalln("get rlimit error: " + err.Error())
	// }
	// rlim.Cur = 1000 * 1000
	// rlim.Max = 1000 * 1000
	// err = syscall.Setrlimit(syscall.RLIMIT_NOFILE, &rlim)
	// if err != nil {
	// 	log.Fatalln("set rlimit error: " + err.Error())
	// }
}

func main() {
	err := RootCmd.Execute()
	if err != nil {
		log.Panic(err)
	}
}
