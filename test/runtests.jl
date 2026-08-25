using Test
using IgBLAST
using Aqua
using CodecZlib

@testset "Aqua.jl" begin
    Aqua.test_all(
        IgBLAST;
        stale_deps=(ignore=[:IgBLAST],),
        deps_compat=(ignore=[:IgBLAST],),
    )
end

@testset "IgBLAST.jl" begin
    @testset "Installation" begin
        @test install_igblast() === nothing
        @test is_igblast_installed() == true
    end

    @testset "native_executable" begin
        @test IgBLAST.native_executable("igblastn"; windows=false) == "igblastn"
        @test IgBLAST.native_executable("igblastn"; windows=true) == "igblastn.exe"
        @test IgBLAST.native_executable("igblastn", Val(false)) == "igblastn"
        @test IgBLAST.native_executable("igblastn", Val(true)) == "igblastn.exe"
    end

    @testset "install_igblast install path (injected)" begin
        fake_sha = Base.SHA1(fill(0x11, 20))
        root = mktempdir()
        bin = joinpath(root, "ncbi-igblast-$(IgBLAST.IGBLAST_VERSION)", "bin")
        mkpath(bin)
        touch(joinpath(bin, IgBLAST.native_executable("igblastn")))

        toml = joinpath(root, "Artifacts.toml")
        write(toml, """
        [IgBLAST]
        git-tree-sha1 = "1111111111111111111111111111111111111111"
        """)

        called = Ref(false)
        fake_ensure = function (name, path)
            called[] = true
            @test name == "IgBLAST"
            @test path == toml
            return nothing
        end

        sha = install_igblast(;
            already_installed=() -> false,
            artifact_toml=toml,
            ensure_fn=fake_ensure,
            artifact_hash_fn=(_, _) -> fake_sha,
            artifact_path_fn=_ -> root,
        )
        @test sha == fake_sha
        @test called[]

        empty_root = mktempdir()
        @test_throws ErrorException install_igblast(;
            already_installed=() -> false,
            artifact_toml=toml,
            ensure_fn=fake_ensure,
            artifact_hash_fn=(_, _) -> fake_sha,
            artifact_path_fn=_ -> empty_root,
        )

        empty_toml = joinpath(empty_root, "empty.toml")
        write(empty_toml, "")
        @test_throws ErrorException install_igblast(;
            already_installed=() -> false,
            artifact_toml=empty_toml,
        )

        rm(root; recursive=true)
        rm(empty_root; recursive=true)
    end

    @testset "verify_igblast_installation" begin
        root = mktempdir()
        sha = Base.SHA1(fill(0x22, 20))
        @test_throws ErrorException IgBLAST.verify_igblast_installation(sha; artifact_path_fn=_ -> root)

        bin = joinpath(root, "ncbi-igblast-$(IgBLAST.IGBLAST_VERSION)", "bin")
        mkpath(bin)
        exe = joinpath(bin, IgBLAST.native_executable("igblastn"))
        touch(exe)
        @test IgBLAST.verify_igblast_installation(sha; artifact_path_fn=_ -> root) == exe
        rm(root; recursive=true)
    end

    @testset "Types and germlines" begin
        @test IgBLASTn <: AbstractIgBLAST
        @test IgBLASTp <: AbstractIgBLAST
        @test IgBLAST.executable(IgBLASTn) == "igblastn"
        @test IgBLAST.executable(IgBLASTp) == "igblastp"
        @test IgBLAST.default_outfmt(IgBLASTn) == 19
        @test IgBLAST.default_outfmt(IgBLASTp) == 7
        @test IgBLAST.supports_auxiliary(IgBLASTn)
        @test !IgBLAST.supports_auxiliary(IgBLASTp)
        @test IgBLAST.auxiliary_capability(IgBLASTn) isa IgBLAST.AuxiliarySupported
        @test IgBLAST.auxiliary_capability(IgBLASTp) isa IgBLAST.AuxiliaryUnsupported
        @test IgBLAST.molecule(IgBLASTn) isa IgBLAST.DNAMolecule
        @test IgBLAST.molecule(IgBLASTp) isa IgBLAST.ProteinMolecule
        @test IgBLAST.blastdb_type(IgBLAST.molecule(IgBLASTn)) == "nucl"
        @test IgBLAST.blastdb_type(IgBLAST.molecule(IgBLASTp)) == "prot"

        @test IgBLAST.normalize_auxiliary(noauxiliary) isa NoAuxiliary
        @test IgBLAST.normalize_auxiliary("aux.txt") isa AuxiliaryFile

        @test germlines_for(IgBLASTn, "V.fa", "D.fa", "J.fa") isa VDJGermlines
        @test germlines_for(IgBLASTp, "V.fa") isa VGermlines
    end

    function write_dummy_fasta(path)
        open(path, "w") do io
            println(io, ">DummySequence")
            println(io, "ACGTACGTACGT")
        end
    end

    function write_protein_query(path)
        open(path, "w") do io
            println(io, ">QuerySequence")
            println(io, "MCRMC")
        end
    end

    function write_nucleotide_db(path)
        open(path, "w") do io
            println(io, ">DummySequence")
            println(io, "ATGCGTATGCGTATGCGT")
        end
    end

    human = Dict{String,String}("organism" => "human", "domain_system" => "imgt")

    @testset "prepare IgBLASTn with custom aux" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        aux_file = tempname() * ".txt"
        output_file = tempname() * ".txt"

        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))
        touch(aux_file)

        result = prepare(
            IgBLASTn, v_database, d_database, j_database;
            aux=aux_file,
            additional_params=Dict{String,String}(
                "organism" => "human",
                "domain_system" => "imgt",
                "ungapped" => "",
            ),
        ) do ig
            @test ig isa IgBLASTSession
            ig(query_file, output_file)
        end

        @test result == output_file
        @test isfile(output_file)
        @test filesize(output_file) > 0

        foreach(rm, (query_file, v_database, d_database, j_database, aux_file, output_file))
    end

    @testset "prepare IgBLASTn without aux" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"

        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))
        dbs = VDJGermlines(v_database, d_database, j_database)

        prepare(IgBLASTn, dbs; additional_params=human) do ig
            @test isconcretetype(typeof(ig))
            ig(query_file, output_file)
        end
        @test isfile(output_file)
        @test filesize(output_file) > 0

        foreach(rm, (query_file, v_database, d_database, j_database, output_file))
    end

    @testset "command is a direct IgBLAST Cmd" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"

        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))
        dbs = VDJGermlines(v_database, d_database, j_database)

        prepare(IgBLASTn, dbs; additional_params=human) do ig
            cmd = @inferred command(ig, query_file, output_file)
            @test cmd isa Cmd
            exe = IgBLAST.executable_path(IgBLASTn)
            @test cmd.exec[1] == exe
            @test "-query" in cmd.exec
            @test "-out" in cmd.exec
            run(cmd)
        end

        @test isfile(output_file)
        @test filesize(output_file) > 0

        foreach(rm, (query_file, v_database, d_database, j_database, output_file))
    end

    @testset "README timing example" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file = tempname() * ".tsv"

        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))
        dbs = VDJGermlines(v_database, d_database, j_database)

        seconds = prepare(IgBLASTn, dbs; additional_params=human, num_threads=1) do ig
            @elapsed run(command(ig, query_file, output_file))
        end

        @test seconds isa Float64
        @test seconds >= 0
        @test isfile(output_file)
        @test filesize(output_file) > 0

        foreach(rm, (query_file, v_database, d_database, j_database, output_file))
    end

    @testset "Missing custom aux errors clearly" begin
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        foreach(write_dummy_fasta, (v_database, d_database, j_database))

        @test_throws ArgumentError prepare(
            IgBLASTn, v_database, d_database, j_database;
            aux="/nonexistent/aux.txt",
        ) do ig
            error("should not run")
        end

        foreach(rm, (v_database, d_database, j_database))
    end

    @testset "IgBLASTp rejects auxiliary data" begin
        v_database = tempname() * ".fasta"
        write_nucleotide_db(v_database)
        aux_file = tempname() * ".txt"
        touch(aux_file)

        @test_throws ArgumentError prepare(IgBLASTp, v_database; aux=aux_file) do ig
            error("should not run")
        end

        rm(v_database)
        rm(aux_file)
    end

    @testset "prepare IgBLASTp (V-only)" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"

        write_protein_query(query_file)
        write_nucleotide_db(v_database)

        prepare(
            IgBLASTp, v_database;
            additional_params=Dict{String,String}("organism" => "human"),
        ) do ig
            ig(query_file, output_file)
        end

        @test isfile(output_file)
        @test filesize(output_file) > 0

        foreach(rm, (query_file, v_database, output_file))
    end

    @testset "Gzip query and output autodetection" begin
        query_file = tempname() * ".fasta"
        query_gz = tempname() * ".fasta.gz"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_gz = tempname() * ".tsv.gz"

        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))
        open(GzipCompressorStream, query_gz, "w") do out
            write(out, read(query_file))
        end

        @test IgBLAST.file_encoding(query_gz) isa IgBLAST.GzipEncoding
        @test IgBLAST.file_encoding(output_gz) isa IgBLAST.GzipEncoding
        @test IgBLAST.file_encoding(query_file) isa IgBLAST.PlainEncoding

        prepare(
            IgBLASTn, v_database, d_database, j_database;
            additional_params=human,
        ) do ig
            ig(query_gz, output_gz)
        end

        @test isfile(output_gz)
        @test filesize(output_gz) > 0
        open(GzipDecompressorStream, output_gz) do io
            text = String(read(io))
            @test !isempty(text)
        end

        foreach(rm, (query_file, query_gz, v_database, d_database, j_database, output_gz))
    end

    @testset "Variant and germline kinds pair by dispatch" begin
        v_database = tempname() * ".fasta"
        write_nucleotide_db(v_database)
        @test_throws MethodError prepare(IgBLASTn, VGermlines(v_database)) do ig
            error("should not run")
        end
        rm(v_database)
    end
end
