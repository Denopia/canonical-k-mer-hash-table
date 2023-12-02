#include "file_reader.hpp"

/*

    File reader class implementation

*/


FastaReader::FastaReader(std::string path, bool reverse) : 
fasta_path(path), current_read_length(0), current_line_number(0), current_read_number(0), reverse_reads_enabled(reverse), current_is_forward(false)
{
    if (!reverse)
    {
        current_is_forward = true;
    }
    current_read = "";
    fasta_file.open(fasta_path);
    if (!fasta_file)
    {
        i_have_reads = false;
    } 
    else
    {
        i_have_reads = true;
        roll_to_next_read();
    }
}

void FastaReader::read_the_next_read()
{
    if (current_line_number == 0)
    {
        std::getline(fasta_file, current_line);
        current_line_number += 1;
    }
    current_read = "";
    while (true) {
        std::getline(fasta_file, current_line);
        current_line_number += 1;
        if (current_line[0] == '>' || current_line[0] == '@' ||current_line.length() == 0)
            break;
        current_read = current_read + current_line;
    }
    current_read_length = current_read.length();
    if (current_read_length > 0)
    {
        current_read_number += 1;
        current_read_length = current_read.length();
    }
    else
    {
        i_have_reads = false;
        fasta_file.close();
    }
}


void FastaReader::roll_to_next_read()
{
    // If we include reverse reads, do this
    if (reverse_reads_enabled){
        // If the previous handled read was a forward read, reverse it
        if (current_is_forward){
            //purestringfunctions::reverse_this_string(current_read);
            current_is_forward = false;
        } else {
            read_the_next_read();
            current_is_forward = true;
        }
    // If we only consider forward reads, simply do this instead
    } else {
        read_the_next_read();
    }
}

int FastaReader::get_current_read_length()
{
    return current_read_length;
}

bool FastaReader::read_is_loaded()
{
    return i_have_reads;
}

char FastaReader::get_current_read_character_at(int position)
{
    if (current_is_forward)
    {
        if (int(current_read_length) > position)
        {
            return current_read.at(position);
        }
        else
        {
            std::cout << "File reader tried to fetch a character when no read was loaded\n";
            exit(1);
        }
    }
    else
    {
        if (int(current_read_length) > position)
        {
            return purestringfunctions::reverse_char(current_read.at(current_read_length - position - 1));
        }
        else
        {
            std::cout << "File reader tried to fetch a character when no read was loaded\n";
            exit(1);
        }
    }
}

std::string FastaReader::get_current_read()
{
    return current_read;
}
