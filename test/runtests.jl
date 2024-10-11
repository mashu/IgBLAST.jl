using Test
using IgBLAST

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
            additional_params = Dict{String, Any}(
                "organism" => "human",
                "domain_system" => "imgt"
            )
        )

        @test isfile(output_file)
        @test filesize(output_file) > 0

        # Clean up temporary files
        for file in [query_file, v_database, d_database, j_database, aux_file, output_file]
            rm(file)
        end
    end
end