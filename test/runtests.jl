using Test
using IgBLAST
using Aqua

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

    @testset "Utility Functions" begin
        @test IgBLAST.get_igblast_url() isa String

        temp_fasta = tempname() * ".fasta"
        open(temp_fasta, "w") do io
            for i in 1:10
                println(io, ">Sequence$i")
                println(io, "ACGT" * repeat("N", i))
            end
        end

        @test IgBLAST.count_fasta_sequences(temp_fasta) == 10
        rm(temp_fasta)
    end

    @testset "Types, traits, and germlines" begin
        @test IgBLASTn <: AbstractIgBLAST
        @test IgBLASTp <: AbstractIgBLAST
        @test IgBLAST.executable(IgBLASTn) == "igblastn"
        @test IgBLAST.executable(IgBLASTp) == "igblastp"
        @test IgBLAST.default_outfmt(IgBLASTn) == 19
        @test IgBLAST.default_outfmt(IgBLASTp) == 7
        @test IgBLAST.supports_auxiliary(IgBLASTn)
        @test !IgBLAST.supports_auxiliary(IgBLASTp)
        @test IgBLAST.molecule(IgBLASTn) isa IgBLAST.DNAMolecule
        @test IgBLAST.molecule(IgBLASTp) isa IgBLAST.ProteinMolecule

        @test IgBLAST.normalize_auxiliary(nothing) isa NoAuxiliary
        @test IgBLAST.normalize_auxiliary("") isa NoAuxiliary
        @test IgBLAST.normalize_auxiliary("aux.txt") isa AuxiliaryFile
        @test IgBLAST.normalize_auxiliary(noauxiliary) isa NoAuxiliary

        @test germlines_for(IgBLASTn, "V.fa", "D.fa", "J.fa") isa VDJGermlines
        @test germlines_for(IgBLASTp, "V.fa", "D.fa", "J.fa") isa VGermlines
        @test germlines_for(IgBLASTp, "V.fa") isa VGermlines
        @test GermlineDatabases === VDJGermlines
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

    @testset "Run IgBLASTn with custom aux" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        aux_file = tempname() * ".txt"
        output_file = tempname() * ".txt"

        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))
        touch(aux_file)

        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTn,
            query_file,
            v_database,
            d_database,
            j_database,
            aux_file,
            output_file;
            additional_params=Dict{String,String}(
                "organism" => "human",
                "domain_system" => "imgt",
                "ungapped" => "",
            ),
        )

        @test isfile(output_file)
        @test filesize(output_file) > 0

        foreach(rm, (query_file, v_database, d_database, j_database, aux_file, output_file))
    end

    @testset "Run IgBLASTn without aux (omitted)" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"

        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))

        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTn,
            query_file,
            v_database,
            d_database,
            j_database,
            output_file;
            additional_params=Dict{String,String}(
                "organism" => "human",
                "domain_system" => "imgt",
            ),
        )

        @test isfile(output_file)
        @test filesize(output_file) > 0

        foreach(rm, (query_file, v_database, d_database, j_database, output_file))
    end

    @testset "Run IgBLASTn without aux (nothing / empty / typed)" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"

        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))
        dbs = VDJGermlines(v_database, d_database, j_database)

        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTn, query_file, dbs, output_file;
            additional_params=Dict{String,String}("organism" => "human", "domain_system" => "imgt"),
        )
        @test isfile(output_file)

        output_file2 = tempname() * ".txt"
        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTn, query_file, v_database, d_database, j_database, nothing, output_file2;
            additional_params=Dict{String,String}("organism" => "human", "domain_system" => "imgt"),
        )
        @test isfile(output_file2)

        output_file3 = tempname() * ".txt"
        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTn, query_file, v_database, d_database, j_database, "", output_file3;
            additional_params=Dict{String,String}("organism" => "human", "domain_system" => "imgt"),
        )
        @test isfile(output_file3)

        foreach(rm, (query_file, v_database, d_database, j_database, output_file, output_file2, output_file3))
    end

    @testset "IgBLASTRunner functor" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"

        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))

        runner = IgBLASTRunner(
            IgBLASTn;
            additional_params=Dict{String,String}(
                "organism" => "human",
                "domain_system" => "imgt",
            ),
        )
        @test runner.aux isa NoAuxiliary
        @test runner.outfmt == 19
        @test_logs (:info, r"IgBLAST analysis completed.*") runner(
            query_file, v_database, d_database, j_database, output_file,
        )
        @test isfile(output_file)

        output_file2 = tempname() * ".txt"
        @test_logs (:info, r"IgBLAST analysis completed.*") runner(
            query_file, VDJGermlines(v_database, d_database, j_database), output_file2,
        )
        @test isfile(output_file2)

        foreach(rm, (query_file, v_database, d_database, j_database, output_file, output_file2))
    end

    @testset "Missing custom aux errors clearly" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"
        foreach(write_dummy_fasta, (query_file, v_database, d_database, j_database))

        @test_throws ArgumentError run_igblast(
            IgBLASTn,
            query_file,
            v_database,
            d_database,
            j_database,
            "/nonexistent/aux.txt",
            output_file,
        )

        foreach(rm, (query_file, v_database, d_database, j_database))
        isfile(output_file) && rm(output_file)
    end

    @testset "Run IgBLASTp (V-only)" begin
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"

        write_protein_query(query_file)
        write_nucleotide_db(v_database)

        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTp,
            query_file,
            v_database,
            output_file;
            additional_params=Dict{String,String}("organism" => "human"),
        )

        @test isfile(output_file)
        @test filesize(output_file) > 0

        # Compat: still accepts unused D/J paths
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file2 = tempname() * ".txt"
        foreach(write_nucleotide_db, (d_database, j_database))
        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTp,
            query_file,
            v_database,
            d_database,
            j_database,
            output_file2;
            additional_params=Dict{String,String}("organism" => "human"),
        )
        @test isfile(output_file2)

        foreach(rm, (query_file, v_database, d_database, j_database, output_file, output_file2))
    end
end
