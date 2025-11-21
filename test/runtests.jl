using Test
using IgBLAST
using Aqua

@testset "Aqua.jl" begin
  Aqua.test_all(
    IgBLAST;
    #ambiguities=(exclude=[SomePackage.some_function], broken=true),
    stale_deps=(ignore=[:IgBLAST],),
    deps_compat=(ignore=[:IgBLAST],),
    #piracies=false,
  )
end
@testset "IgBLAST.jl" begin
    @testset "Installation" begin
        @test install_igblast() === nothing
        @test is_igblast_installed() == true
    end

    @testset "Utility Functions" begin
        @test IgBLAST.get_igblast_url() isa String

        # Create a temporary FASTA file for testing
        temp_fasta = tempname() * ".fasta"
        open(temp_fasta, "w") do io
            for i in 1:10
                println(io, ">Sequence$i")
                println(io, "ACGT" * repeat("N", i))
            end
        end

        @test IgBLAST.count_fasta_sequences(temp_fasta) == 10

        # Clean up the temporary file
        rm(temp_fasta)
    end

    @testset "IgBLAST Types" begin
        @test IgBLASTn <: AbstractIgBLAST
        @test IgBLASTp <: AbstractIgBLAST
        @test IgBLAST.executable(IgBLASTn) == "igblastn"
        @test IgBLAST.executable(IgBLASTp) == "igblastp"
    end

    @testset "Run IgBLAST" begin
        # Create temporary files for testing
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        aux_file = tempname() * ".txt"
        output_file = tempname() * ".txt"
        # Write some dummy content to the files
        for file in [query_file, v_database, d_database, j_database]
            open(file, "w") do io
                println(io, ">DummySequence")
                println(io, "ACGTACGTACGT")
            end
        end
        touch(aux_file)

        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTn,
            query_file,
            v_database,
            d_database,
            j_database,
            aux_file,
            output_file,
            additional_params = Dict{String, String}(
                "organism" => "human",
                "domain_system" => "imgt",
                "ungapped" => ""
                )
        )

        @test isfile(output_file)
        @test filesize(output_file) > 0

        # Clean up temporary files
        for file in [query_file, v_database, d_database, j_database, aux_file, output_file]
            rm(file)
        end
    end

    @testset "Run IgBLASTn without aux_file" begin
        # Test IgBLASTn with empty aux_file (optional for assignments without CDR3)
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"
        
        # Write some dummy content to the files
        for file in [query_file, v_database, d_database, j_database]
            open(file, "w") do io
                println(io, ">DummySequence")
                println(io, "ACGTACGTACGT")
            end
        end

        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTn,
            query_file,
            v_database,
            d_database,
            j_database,
            "",  # Empty aux_file
            output_file,
            additional_params = Dict{String, String}(
                "organism" => "human",
                "domain_system" => "imgt"
            )
        )

        @test isfile(output_file)
        @test filesize(output_file) > 0

        # Clean up temporary files
        for file in [query_file, v_database, d_database, j_database, output_file]
            rm(file)
        end
    end

    @testset "Run IgBLASTp" begin
        # Create temporary files for testing IgBLASTp
        # For IgBLASTp: query should be protein sequences, database should be nucleotide sequences
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        aux_file = tempname() * ".txt"
        output_file = tempname() * ".txt"
        
        # Query file: protein sequences
        open(query_file, "w") do io
            println(io, ">QuerySequence")
            println(io, "MCRMC")  # Protein sequence
        end
        
        # Database files: nucleotide sequences (will be translated to protein during DB preparation)
        # Using sequences that are multiples of 3 for proper translation
        for file in [v_database, d_database, j_database]
            open(file, "w") do io
                println(io, ">DummySequence")
                println(io, "ATGCGTATGCGTATGCGT")  # 18 nucleotides = 6 amino acids (MCRMCR)
            end
        end
        touch(aux_file)

        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTp,
            query_file,
            v_database,
            d_database,
            j_database,
            aux_file,
            output_file,
            additional_params = Dict{String, String}(
                "organism" => "human"
            )
        )

        @test isfile(output_file)
        @test filesize(output_file) > 0

        # Clean up temporary files
        for file in [query_file, v_database, d_database, j_database, aux_file, output_file]
            rm(file)
        end
    end

    @testset "Run IgBLASTp without aux_file" begin
        # Test IgBLASTp with empty aux_file
        query_file = tempname() * ".fasta"
        v_database = tempname() * ".fasta"
        d_database = tempname() * ".fasta"
        j_database = tempname() * ".fasta"
        output_file = tempname() * ".txt"
        
        # Query file: protein sequences
        open(query_file, "w") do io
            println(io, ">QuerySequence")
            println(io, "MCRMC")  # Protein sequence
        end
        
        # Database files: nucleotide sequences (will be translated to protein during DB preparation)
        for file in [v_database, d_database, j_database]
            open(file, "w") do io
                println(io, ">DummySequence")
                println(io, "ATGCGTATGCGTATGCGT")  # 18 nucleotides = 6 amino acids
            end
        end

        @test_logs (:info, r"IgBLAST analysis completed.*") run_igblast(
            IgBLASTp,
            query_file,
            v_database,
            d_database,
            j_database,
            "",  # Empty aux_file
            output_file,
            additional_params = Dict{String, String}(
                "organism" => "human"
            )
        )

        @test isfile(output_file)
        @test filesize(output_file) > 0

        # Clean up temporary files
        for file in [query_file, v_database, d_database, j_database, output_file]
            rm(file)
        end
    end
end
