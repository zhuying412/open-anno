package prepare

import "github.com/spf13/cobra"

func NewPrepareCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "prepare",
		Short: "Prepare required database files",
	}
	cmd.AddCommand(NewPrepareTranscriptCmd())
	cmd.AddCommand(NewPrepareDatabaseCmd())
	return cmd
}
