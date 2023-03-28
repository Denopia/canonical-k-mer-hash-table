#include "program_runs.hpp"

/*
    Run canonical k-mer hash table mode 
*/
int run_mode_1(int argc, char const* argv[])
{
    //uint64_t b0 = 31/32;
    //uint64_t b1 = 32/32;
    //uint64_t b2 = 33/32;
    //std::cout << b0 << " " <<  b1 << " " << b2 << "\n"; 
    //std::cout << "Canonical hash table mode run starts\n";

    // --- Initialize needed arguments ---

    // k-mer length
    int k = 0;
    // hash table slots
    int min_slots = 0;
    // max number of uniq entries in the hash table
    int uniq_entries_limit = 0;
    // run queries after hash map is built
    bool query_enabled = false;
    // first full k-mer of the read is checked before anything else
    bool first_kmer_check_enabled = false;
    // reads file path
    std::string reads_path = "No file";
    // Help trigger
    bool user_wants_help = false;
    // Verbosity flag
    bool verbose = false;

    int min_abundance = 1;
    bool write_output = false;
    std::string output_path = "nessu";
    std::string output_directory = ".";

    int argi = 1;

    // --- Parse arguments ---
    while (argi < argc)
    {
        std::string as(argv[argi]);
        
        if (as.compare("-k") == 0){
            k = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-s") == 0)
        {
            min_slots = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-u") == 0)
        {
            uniq_entries_limit = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-q") == 0)
        {
            query_enabled = true;
            argi += 1;
        }
        else if (as.compare("-f") == 0)
        {
            first_kmer_check_enabled = true;
            argi += 1;
        }
        else if (as.compare("-p") == 0)
        {
            reads_path = std::string(argv[argi+1]);
            argi += 2;
        }
        else if (as.compare("-h") == 0)
        {
            user_wants_help = true;
            argi += 1;
        }
        else if (as.compare("-v") == 0)
        {
            verbose = true;
            argi += 1;
        }
        else if (as.compare("-a") == 0){
            min_abundance = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-o") == 0)
        {
            output_path = std::string(argv[argi+1]);
            write_output = true;
            argi += 2;
        }
        else if (as.compare("-d") == 0)
        {
            output_directory = std::string(argv[argi+1]);
            argi += 2;
        }
    }

    if (user_wants_help)
    {   std::cout << "This is a k-mer has table program. Here are the accepted arguments:\n\n";
        std::cout << "\t-k [integer] determines the k-mer length\n";
        std::cout << "\t-a [integer] minimum k-mer abundance (for writing output, default: 1)\n";
        std::cout << "\t-s [integer] determines the minimum hash table size (actual used size is the next prime)\n";
        std::cout << "\t-u [integer] determines a limit for the number of items that can be stored in the hash table (default: maximum size)\n";
        std::cout << "\t-q flag enables queries after hash table is built\n";
        std::cout << "\t-f flag tells the program to do a first k-mer check after buffer is empty (should reduce the number of stored entries)\n";
        std::cout << "\t-p [string] indicates the path to the read file (at the moment fasta format is required)\n";
        std::cout << "\t-v flag tells the program to print what is going on during the program\n";
        std::cout << "\t-o [string] tells the program to write found k-mers into an output file with this name\n";
        std::cout << "\t-d [string] shows the directory where output is written\n";
        std::cout << "\t-h flag shows this help message\n\n";

        return 0;
    }

    // --- Next, verify some arguments ---    
     if (reads_path == "No File")
    {
        std::cout << "** ERROR ** I need a path to the reads file (must be fasta)...\n";
        return 1;
    }
    if (k < 1)
    {
        std::cout << "** ERROR ** The value of k must be greater than 0...\n";
        return 1;
    }
    if (min_slots == 0)
    {
        std::cout << "** ERROR ** Minimum number of hash table slots not given\n";
        return 1;
        
    }

    // --- Make the number of hash table slots prime ---
    uint32_t prime_slots = mathfunctions::next_prime(min_slots);
    
    // --- Define new variable for hash table size ---
    int hash_table_slots = prime_slots;

    // --- Set the desired load factor
    float load_factor = 0.8;

    // --- Build file reader ---
    bool reverse_reads_enabled = false;
    FastaReader* file_reader = new FastaReader(reads_path, reverse_reads_enabled);

    // --- Build k-mer factory ---
    KMerFactoryCanonical2BC* kmer_factory = new KMerFactoryCanonical2BC(k);
    
    // --- Build rolling hasher ----
    RollingHasherDual* hasher = new RollingHasherDual(hash_table_slots, k);

    // --- Build hash table ---
    PointerHashTableCanonical* hash_table = new PointerHashTableCanonical(hash_table_slots, k, kmer_factory->number_of_blocks);

    uint64_t current_kmer_slot = 0;
    uint64_t predecessor_kmer_slot = hash_table_slots;

    uint64_t current_kmer_hash = 0;

    bool predecessor_kmer_exists = false;

    int read_position = 0;
    int first_check_read_position = 0;
    bool first_kmer_unread = true;

    int reads_processed = 0;

    // --- Counters for new k-mer statistics ---
    uint64_t counter_10k = 0;
    uint64_t counter_100k = 0;
    uint64_t counter_1m = 0;
    uint64_t counter_10m = 0;

    uint64_t added_kmers_at_last_10k_mark = 0;
    uint64_t added_kmers_at_last_100k_mark = 0;
    uint64_t added_kmers_at_last_1m_mark = 0;
    uint64_t added_kmers_at_last_10m_mark = 0;

    uint64_t last_10k_interval_new_kmers = 0;
    uint64_t last_100k_interval_new_kmers = 0;
    uint64_t last_1m_interval_new_kmers = 0;
    uint64_t last_10m_interval_new_kmers = 0;


    // --- Starting reading the reads and adding found k-mers to the hash table ---
    while (file_reader->read_is_loaded())
    {   
        // if read is too short it is skipped
        if (file_reader->get_current_read_length() < k)
        {
            file_reader->roll_to_next_read();
            reads_processed += 1;
            continue;
        }

        predecessor_kmer_exists = false;
        predecessor_kmer_slot = hash_table_slots;
        read_position = 0;

        //std::cout << "=======================================================================\n";
        //std::cout << "============ This is read number " << reads_processed << " ============\n";
        //std::cout << "=======================================================================\n";

        // Start reading characters from the read until it is fully read
        while (read_position < file_reader->get_current_read_length())
        {
            // Push character into the k-mer factory buffer
            kmer_factory->push_new_character(file_reader->get_current_read_character_at(read_position));
            

            // If it was not a legit character, reset hahser too
            if (kmer_factory->get_number_of_stored_characters() == 0){
                hasher->reset();
                predecessor_kmer_exists = false;
                predecessor_kmer_slot = hash_table_slots;
            } else {
                hasher->update_rolling_hash(kmer_factory->get_forward_newest_character(), kmer_factory->get_forward_pushed_off_character());
            }

            //std::cout << "---------- Read position is " << read_position << " ----------\n";
            //std::cout << "File reader character is " << file_reader->get_current_read_character_at(read_position) << "\n";
            //std::cout << "Newest character is "  << kmer_factory->get_forward_newest_character() << "\n";
            //std::cout << "Oldest character is "  << kmer_factory->get_forward_pushed_off_character() << "\n";
            //std::cout << "Forward factory block is: " << kmer_factory->get_forward_block(0) << "\n";
            //std::cout << "Backward factory block is: " << kmer_factory->get_backward_block(0) << "\n";
            //std::cout << "k-mer hash forward is " << hasher->get_current_hash_forward() << "\n";
            //std::cout << "k-mer hash backward is " << hasher->get_current_hash_backward() << "\n\n";

            // If k-mer factory has k characters stored, update the hash table
            if (kmer_factory->get_number_of_stored_characters() == k)
            //if ((kmer_factory->get_number_of_stored_characters() == k) || (true))
            {
                // Add the current k-mer in the k-mer factory buffer
                // If it is already in the buffer, count is increased by one
                // If it is not, k-mer is inserted in the hash table with count 1
                current_kmer_slot = hash_table->process_kmer(kmer_factory, hasher, predecessor_kmer_exists, predecessor_kmer_slot); 
                predecessor_kmer_exists = true;
                predecessor_kmer_slot = current_kmer_slot;

                // Update counters to track how many new k-mers have ben added in the past
                counter_10k += 1;
                counter_100k += 1;
                counter_1m += 1;
                counter_10m += 1;

                /*
                if (counter_10k == 10000)
                {
                    last_10k_interval_new_kmers = hash_table->get_number_of_inserted_items() - added_kmers_at_last_10k_mark;
                    added_kmers_at_last_10k_mark = hash_table->get_number_of_inserted_items();
                    counter_10k = 0;
                    if (counter_100k == 100000)
                    {
                        last_100k_interval_new_kmers = hash_table->get_number_of_inserted_items() - added_kmers_at_last_100k_mark;
                        added_kmers_at_last_100k_mark = hash_table->get_number_of_inserted_items();
                        counter_100k = 0;
                        if (counter_1m == 1000000)
                        {
                            last_1m_interval_new_kmers = hash_table->get_number_of_inserted_items() - added_kmers_at_last_1m_mark;
                            added_kmers_at_last_1m_mark = hash_table->get_number_of_inserted_items();
                            counter_1m = 0;
                            if (counter_10m == 10000000)
                            {
                                last_10m_interval_new_kmers = hash_table->get_number_of_inserted_items() - added_kmers_at_last_10m_mark;
                                added_kmers_at_last_10m_mark = hash_table->get_number_of_inserted_items();
                                counter_10m = 0;
                            }
                        }
                    }
                }
                */
            }

            // If hash table is too full, kill program
            
            if (hash_table->get_number_of_inserted_items() > load_factor*hash_table_slots)
            {
                std::cout << "Number of k-mers exceeds desired load factor. Aborting program.\n";
                return 1;
            }
            
            read_position+=1;

        }

        file_reader->roll_to_next_read();
        kmer_factory->reset();
        hasher->reset();

        reads_processed += 1;

        
        if (verbose && (reads_processed % 100 == 0))
        {
            std::cout << "Processed read number " << reads_processed << "\n";
            //std::cout << hash_table->get_number_of_inserted_complete_kmers() << "/" << hash_table_slots << " k-mers stored currently in main array\n";
            std::cout << hash_table->get_number_of_inserted_items() << "/" << hash_table_slots << " k-mers stored currently\n";
            std::cout << hash_table->get_number_of_inserted_items_in_main() << "/" << hash_table_slots << " k-mers stored currently in main array\n";
            std::cout << "--------------------------------------------------------\n";
            std::cout << hash_table->get_max_number_of_secondary_slots_in_use() << "/" << hash_table->get_number_of_max_secondary_slots() << " slots touched in secondary array\n";
            std::cout << hash_table->get_number_of_secondary_slots_in_use() << "/" << hash_table->get_number_of_max_secondary_slots() << " k-mers stored currently in secondary array\n";
            //std::cout << "--------------------------------------------------------\n";
            //std::cout << "This many new k-mers were added in the latest 10k k-mer interval:  " << last_10k_interval_new_kmers / 10000.0 << "\n";
            //std::cout << "This many new k-mers were added in the latest 100k k-mer interval: " << last_100k_interval_new_kmers / 100000.0 << "\n";
            //std::cout << "This many new k-mers were added in the latest 1m k-mer interval:   " << last_1m_interval_new_kmers / 1000000.0 << "\n";
            //std::cout << "This many new k-mers were added in the latest 10m k-mer interval:  " << last_10m_interval_new_kmers / 10000000.0 << "\n";
            std::cout << "===============================================================================\n";
        }

    }

    
    if (write_output)
    {
        uint64_t writing_microseconds;
        std::chrono::high_resolution_clock::time_point start_time, end_time;
        
        /*
        start_time = std::chrono::high_resolution_clock::now();

        std::cout << "Writing k-mers to output file\n";
        hash_table->write_kmers_on_disk_separately(min_abundance, output_path);

        end_time = std::chrono::high_resolution_clock::now();

        writing_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time).count();

        std::cout << "### Writing time: " << writing_microseconds << " microseconds\n";


        start_time = std::chrono::high_resolution_clock::now();

        std::cout << "Writing k-mers to output file faster\n";
        hash_table->write_kmers_on_disk_separately_faster(min_abundance, output_path+"FASTER.txt");

        end_time = std::chrono::high_resolution_clock::now();

        writing_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time).count();

        std::cout << "### Writing time faster: " << writing_microseconds << " microseconds\n";
        
        start_time = std::chrono::high_resolution_clock::now();

        std::cout << "Writing k-mers to output file even faster\n";
        hash_table->write_kmers_on_disk_separately_even_faster(min_abundance, output_path+"EVENFASTER.txt");

        end_time = std::chrono::high_resolution_clock::now();

        writing_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time).count();

        std::cout << "### Writing time even faster: " << writing_microseconds << " microseconds\n";
        */
        start_time = std::chrono::high_resolution_clock::now();

        std::cout << "Writing k-mers to output file even faster\n";
        hash_table->write_kmers_on_disk_separately_even_faster(min_abundance, output_path);

        end_time = std::chrono::high_resolution_clock::now();

        writing_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time).count();

        std::cout << "### Writing time even faster: " << writing_microseconds << " microseconds\n";
        
    }


    // --- Run queries if they are enabled ---
    
    if (query_enabled)
    {
        std::string query_kmer;
        if (verbose)
            std::cout << "Query k-mers (input q to quit):\n";
        while (true)
        {
            current_kmer_slot = 0;

            std::cin >> query_kmer;
            if (!query_kmer.compare("q") || !query_kmer.compare("Q"))
            {
                if (verbose)
                    std::cout << "Query end now.\n";
                break;
            }    
            else if (query_kmer.length()!=uint64_t(k))
                if (verbose)
                    std::cout << "k-mer must be " << k << " long...\n";
                else
                    std::cout << "-1\n";
            else if (!twobitstringfunctions::is_clean_string(query_kmer))
                if (verbose)
                    std::cout << "Your query k-mer is not valid...\n";
                else
                    std::cout << "-1\n";
            else 
            {
                kmer_factory->reset();
                hasher->reset();

                for (char cc : query_kmer)
                {
                    kmer_factory->push_new_character(cc);
                    hasher->update_rolling_hash(kmer_factory->get_forward_newest_character(), kmer_factory->get_forward_pushed_off_character());
                }
                current_kmer_slot = hash_table->find(kmer_factory, hasher, false, hash_table_slots);
                if (kmer_factory->forward_kmer_is_canonical())
                {
                    std::cout << "Querying for " << query_kmer << "\n";
                }
                else 
                {
                    std::string qkr = query_kmer;
                    purestringfunctions::reverse_this_string(qkr);
                    std::cout << "Querying for " << qkr << " (the canonical form of " << query_kmer << ")\n";
                }
                if (current_kmer_slot == hash_table_slots){
                    std::cout << "0\n";
                } else {
                    std::cout << hash_table->get_kmer_count_in_slot(current_kmer_slot) << "\n";
                }
            }
        }
    }
    

    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int countM = 0;
    
    for (uint32_t l = 0; l < hash_table_slots; l++){
        if (hash_table->slot_is_occupied(l)){
            if (hash_table->get_kmer_count_in_slot(l) == 1){
                count1 += 1;
            } else if (hash_table->get_kmer_count_in_slot(l) == 2){
                count2 += 1;
            } else if (hash_table->get_kmer_count_in_slot(l) == 3){
                count3 += 1;
            } else if (hash_table->get_kmer_count_in_slot(l) > 3){
                countM += 1;
            }
        }
    }
    

    std::cout << "- ABUNDANCE COUNTS -\n";
    std::cout << "k-mer abundance = 1 : " << count1 << "\n";
    std::cout << "k-mer abundance = 2 : " << count2 << "\n";
    std::cout << "k-mer abundance = 3 : " << count3 << "\n";
    std::cout << "k-mer abundance > 3 : " << countM << "\n";

    std::cout << "- CUMULATIVE ABUNDANCE COUNTS -\n";
    std::cout << "k-mer abundance >= 1 : " << count1+count2+count3+countM << "\n";
    std::cout << "k-mer abundance >= 2 : " << count2+count3+countM << "\n";
    std::cout << "k-mer abundance >= 3 : " << count3+countM << "\n";
    std::cout << "k-mer abundance >  3 : " << countM << "\n";

    // --- Free allocated memory ---
    //delete file_reader;
    //delete kmer_factory;
    //delete hasher;
    // Hash table was not deleted in previous version ?? 
    //delete hash_table; 

    std::cout << "Number of stored k-mers:  " << hash_table->get_number_of_inserted_items() << "\n";
    // --- Program run is now finished ---
    std::cout << "Canonical hash table mode run ends\n";

    return 0;

}
