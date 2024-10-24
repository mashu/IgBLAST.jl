var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = IgBLAST","category":"page"},{"location":"#IgBLAST","page":"Home","title":"IgBLAST","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for IgBLAST.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [IgBLAST]","category":"page"},{"location":"#IgBLAST.IgBLAST","page":"Home","title":"IgBLAST.IgBLAST","text":"IgBLAST\n\nA Julia package for running IgBLAST analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.\n\nThis package provides a convenient interface to install and run IgBLAST, supporting both IgBLASTn and IgBLASTp variants. It handles the installation of IgBLAST binaries, prepares input files, runs analyses, and monitors progress.\n\nExports\n\ninstall_igblast: Function to install IgBLAST binaries\nrun_igblast: Main function to run IgBLAST analyses\nis_igblast_installed: Function to check if IgBLAST is installed\nIgBLASTn: Type representing the nucleotide version of IgBLAST\nIgBLASTp: Type representing the protein version of IgBLAST\n\nExamples\n\nusing IgBLAST\n\n# Install IgBLAST if not already installed\ninstall_igblast()\n\n# Run an IgBLASTn analysis\nrun_igblast(\n    IgBLASTn,\n    \"query.fasta\",\n    \"V.fasta\",\n    \"D.fasta\",\n    \"J.fasta\",\n    \"auxiliary.txt\",\n    \"output.txt\",\n    additional_params = Dict(\"organism\" => \"human\", \"domain_system\" => \"imgt\")\n)\n\n\n\n\n\n","category":"module"},{"location":"#IgBLAST.AbstractIgBLAST","page":"Home","title":"IgBLAST.AbstractIgBLAST","text":"AbstractIgBLAST\n\nAn abstract type representing IgBLAST variants.\n\n\n\n\n\n","category":"type"},{"location":"#IgBLAST.IgBLASTn","page":"Home","title":"IgBLAST.IgBLASTn","text":"IgBLASTn\n\nA concrete type representing the nucleotide version of IgBLAST.\n\n\n\n\n\n","category":"type"},{"location":"#IgBLAST.IgBLASTp","page":"Home","title":"IgBLAST.IgBLASTp","text":"IgBLASTp\n\nA concrete type representing the protein version of IgBLAST.\n\n\n\n\n\n","category":"type"},{"location":"#IgBLAST.executable-Tuple{Type{IgBLASTn}}","page":"Home","title":"IgBLAST.executable","text":"executable(::Type{T}) where T <: AbstractIgBLAST\n\nReturns the executable name for a given IgBLAST variant.\n\n\n\n\n\n","category":"method"},{"location":"#IgBLAST.run_igblast-Union{Tuple{T}, Tuple{Type{T}, Vararg{String, 6}}, Tuple{Type{T}, String, String, String, String, String, String, Int64}} where T<:AbstractIgBLAST","page":"Home","title":"IgBLAST.run_igblast","text":"run_igblast(\n    igblast_type::Type{T},\n    query_file::String,\n    v_database::String,\n    d_database::String,\n    j_database::String,\n    aux_file::String,\n    output_file::String,\n    num_threads::Int = Base.Threads.nthreads();\n    outfmt::Int = 19,\n    additional_params::Dict{String, Any} = Dict()\n) where T <: AbstractIgBLAST\n\nRun IgBLAST with the specified parameters.\n\nRequired parameters:\n\nigblast_type: Type of IgBLAST to run (IgBLASTn or IgBLASTp)\nquery_file: Path to the input query file\nv_database: Path to the V gene database\nd_database: Path to the D gene database\nj_database: Path to the J gene database\naux_file: Path to the auxiliary file\noutput_file: Path for the output file\n\nOptional parameters:\n\nnum_threads: Number of threads to use (default: all available threads)\noutfmt: Output format (default: 19 for AIRR format)\nadditional_params: Dictionary of additional IgBLAST parameters\n\nExample:\n\nrun_igblast(IgBLASTn, \"query.fasta\", \"V.fasta\", \"D.fasta\", \"J.fasta\", \"aux.txt\", \"output.txt\",\n            additional_params = Dict(\n                \"organism\" => \"human\",\n                \"domain_system\" => \"imgt\"\n            ))\n\n\n\n\n\n","category":"method"}]
}
