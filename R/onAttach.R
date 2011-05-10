.onAttach <- function(lib, pkg) {
	meta <- packageDescription("mvmeta")
	attachmsg <- paste("Package mvmeta ",meta$Version,
		". For an overview type: help('mvmeta-package').",sep="")
	packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
}