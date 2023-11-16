/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#if 0
#include <cmath>
#include "sav/import.hpp"
#include "sav/sort.hpp"
#include "sav/utility.hpp"
#include "savvy/reader.hpp"
#include "savvy/savvy.hpp"
#include "savvy/writer.hpp"

#include <cstdlib>
#include <getopt.h>

#include <fstream>
#include <vector>
#include <set>

class import_prog_args
{
private:
  std::vector<option> long_options_;
  std::unordered_set<std::string> subset_ids_;
  std::vector<savvy::genomic_region> regions_;
  std::unordered_set<std::string> pbwt_fields_;
  std::unordered_set<std::string> sparse_fields_ = {"GT", "HDS", "EC", "DS"};
  std::string input_path_;
  std::string output_path_;
  std::string index_path_;
  savvy::phasing phasing_ = savvy::phasing::unknown;
  double sparse_threshold_ = 1.0;
  int update_info_ = -1;
  int compression_level_ = -1;
  std::uint16_t block_size_ = savvy::writer::default_block_size;
  bool help_ = false;
  bool index_ = false;
  savvy::fmt format_ = savvy::fmt::gt;
  savvy::bounding_point bounding_point_ = savvy::bounding_point::beg;
  std::unique_ptr<savvy::s1r::sort_point> sort_type_ = nullptr;
  savvy::vcf::empty_vector_policy empty_vector_policy_ = savvy::vcf::empty_vector_policy::fail;
public:
  import_prog_args() :
    long_options_(
      {
        {"block-size", required_argument, 0, 'b'},
        {"bounding-point", required_argument, 0, 'p'},
        {"data-format", required_argument, 0, 'd'},
        {"help", no_argument, 0, 'h'},
        {"index", no_argument, 0, 'x'},
        {"index-file", required_argument, 0, 'X'},
        {"phasing", required_argument, 0, '\x01'},
        {"pbwt-fields", required_argument, 0, '\x01'},
        {"regions", required_argument, 0, 'r'},
        {"regions-file", required_argument, 0, 'R'},
        {"sample-ids", required_argument, 0, 'i'},
        {"sample-ids-file", required_argument, 0, 'I'},
        {"skip-empty-vectors", no_argument, 0, '\x01'},
        {"sort", no_argument, 0, 's'},
        {"sort-point", required_argument, 0, 'S'},
        {"sparse-fields", required_argument, 0, '\x01'},
        {"sparse-threshold", required_argument, 0, '\x01'},
        {"update-info", required_argument, 0, '\x01'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  const std::string& index_path() const { return index_path_; }
  const std::unordered_set<std::string>& subset_ids() const { return subset_ids_; }
  const std::unordered_set<std::string>& pbwt_fields() const { return pbwt_fields_; }
  const std::unordered_set<std::string>& sparse_fields() const { return sparse_fields_; }
  const std::vector<savvy::genomic_region>& regions() const { return regions_; }
  std::uint8_t compression_level() const { return std::uint8_t(compression_level_); }
  std::uint16_t block_size() const { return block_size_; }
  savvy::fmt format() const { return format_; }
  savvy::bounding_point bounding_point() const { return bounding_point_; }
  const std::unique_ptr<savvy::s1r::sort_point>& sort_type() const { return sort_type_; }
  savvy::vcf::empty_vector_policy empty_vector_policy() const { return empty_vector_policy_; }
  savvy::phasing phasing() { return phasing_; }
  double sparse_threshold() { return sparse_threshold_; }
  bool update_info() const { return update_info_ != 0; }
  bool index_is_set() const { return index_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os) const
  {
    os << "Usage: sav import [opts ...] [in.{vcf,vcf.gz,bcf}] [out.sav]\n";
    os << "\n";
    os << " -#                        Number (#) of compression level (1-19; default: " << savvy::writer::default_compression_level << ")\n";
    os << " -b, --block-size          Number of markers in compression block (0-65535; default: " << savvy::writer::default_block_size << ")\n";
    os << " -d, --data-format         Format field to copy (GT or HDS, default: GT)\n";
    os << " -h, --help                Print usage\n";
    os << " -i, --sample-ids          Comma separated list of sample IDs to subset\n";
    os << " -I, --sample-ids-file     Path to file containing list of sample IDs to subset\n";
    os << " -p, --bounding-point      Determines the inclusion policy of indels during region queries (any, all, beg, or end; default is beg)\n";
    os << " -r, --regions             Comma separated list of regions formated as chr[:start-end]\n";
    os << " -R, --regions-file        Path to file containing list of regions formatted as chr<tab>start<tab>end\n";
    os << " -s, --sort                Enables sorting by first position of allele\n";
    os << " -S, --sort-point          Enables sorting and specifies which allele position to sort by (beg, mid, or end)\n";
    os << " -x, --index               Enables indexing\n";
    os << " -X, --index-file          Enables indexing and specifies index output file\n";
    os << "\n";
    os << "     --phasing             Sets file phasing status if phasing header is not present (none, full, or partial)\n";
    os << "     --pbwt-fields         Comma separated list of FORMAT fields for which to enable PBWT sorting\n";
    os << "     --skip-empty-vectors  Skips variants that don't contain the request data format (By default, the import fails)\n";
    os << "     --sparse-fields       Comma separated list of FORMAT fields to make sparse (default: GT,HDS,DS,EC)\n";
    os << "     --sparse-threshold    Non-zero frequency threshold for which sparse fields are encoded as sparse vectors (default: 1.0)\n";
    os << "     --update-info         Specifies whether AC, AN, AF and MAF info fields should be updated (always, never or auto, default: auto)\n";

    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "0123456789b:d:f:hi:I:p:r:R:sS:xX:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
        case '\x01':
        {
          if (strcmp(long_options_[long_index].name, "skip-empty-vectors") == 0)
          {
            empty_vector_policy_ = savvy::vcf::empty_vector_policy::skip;
            break;
          }
          else if (strcmp(long_options_[long_index].name, "phasing") == 0)
          {
            if (strcmp(optarg, "full") == 0)
              phasing_ = savvy::phasing::full;
            else if (strcmp(optarg, "none") == 0)
              phasing_ = savvy::phasing::none;
            else if (strcmp(optarg, "partial") == 0)
              phasing_ = savvy::phasing::partial;
            else
            {
              std::cerr << "Invalid --phasing argument (" << optarg << ")\n";
              return false;
            }
            break;
          }
          else if (strcmp(long_options_[long_index].name, "pbwt-fields") == 0)
          {
            pbwt_fields_ = split_string_to_set(optarg, ',');
            break;
          }
          else if (strcmp(long_options_[long_index].name, "sparse-fields") == 0)
          {
            sparse_fields_ = split_string_to_set(optarg, ',');
            break;
          }
          else if (strcmp(long_options_[long_index].name, "sparse-threshold") == 0)
          {
            sparse_threshold_ = std::atof(optarg);
            break;
          }
          std::cerr << "Invalid long only index (" << long_index << ")\n";
          return false;
        }
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          if (compression_level_ < 0)
            compression_level_ = 0;
          compression_level_ *= 10;
          compression_level_ += copt - '0';
          break;
        case 'b':
          block_size_ = std::uint16_t(std::atoi(optarg) > 0xFFFF ? 0xFFFF : std::atoi(optarg));
          break;
        case 'd':
        {
          std::string str_opt_arg(optarg ? optarg : "");
          if (str_opt_arg == "HDS")
          {
            format_ = savvy::fmt::hds;
          }
          else if (str_opt_arg != "GT")
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        }
        case 'h':
          help_ = true;
          return true;
        case 'r':
          for (const auto& r : split_string_to_vector(optarg, ','))
            regions_.emplace_back(string_to_region(r));
          break;
        case 'R':
          for (const auto& r : split_file_to_vector(optarg))
          {
            std::string s = r;
            std::size_t pos = s.find('\t');
            if (pos != std::string::npos)
            {
              s[pos] = ':';
              pos = s.find('\t', pos + 1);
              if (pos != std::string::npos)
                s[pos] = '-';
            }
            regions_.emplace_back(string_to_region(s));
          }
          break;
        case 'i':
          subset_ids_ = split_string_to_set(optarg, ',');
          break;
        case 'I':
          subset_ids_ = split_file_to_set(optarg);
          break;
        case 'p':
        {
          std::string str_opt_arg(optarg ? optarg : "");
          if (str_opt_arg == "any")
          {
            bounding_point_ = savvy::bounding_point::any;
          }
          else if (str_opt_arg == "all")
          {
            bounding_point_ = savvy::bounding_point::all;
          }
          else if (str_opt_arg == "beg")
          {
            bounding_point_ = savvy::bounding_point::beg;
          }
          else if (str_opt_arg == "end")
          {
            bounding_point_ = savvy::bounding_point::end;
          }
          else
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        }
        case 's':
          sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::beg);
          break;
        case 'S':
        {
          std::string sort_str(optarg);
          if (sort_str.size())
          {
            if (sort_str.front()=='b')
            {
              sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::beg);
            }
            else if (sort_str.front()=='e')
            {
              sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::end);
            }
            else if (sort_str.front()=='m')
            {
              sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::mid);
            }
            else
            {
              std::cerr << "Invalid --sort-point argument (" << sort_str << ")." << std::endl;
              return false;
            }
          }
          break;
        }
        case 'x':
          index_ = true;
          break;
        case 'X':
          index_ = true;
          index_path_ = optarg;
          break;
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count < 2 && index_ && index_path_.empty())
    {
      std::cerr << "--index-file must be specified when output path is not." << std::endl;
      return false;
    }

    if (remaining_arg_count == 0)
    {
      if (regions_.size())
      {
        std::cerr << "Input path must be specified when using --regions option." << std::endl;
        return false;
      }

      input_path_ = "/dev/stdin";
      output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
      output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 2)
    {
      input_path_ = argv[optind];
      output_path_ = argv[optind + 1];

      if (index_ && index_path_.empty())
        index_path_ = output_path_ + ".s1r";
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    if (update_info_ < 0)
    {
      update_info_ = subset_ids_.size() ? 1 : 0; // Automatically update info fields if samples are subset.
    }

    if (compression_level_ < 0)
      compression_level_ = savvy::writer::default_compression_level;
    else if (compression_level_ > 19)
      compression_level_ = 19;

    return true;
  }
};


int import_records(savvy::vcf::indexed_reader<1>& in, const std::vector<savvy::genomic_region>& regions, savvy::fmt data_format, bool update_info, savvy::sav::writer& out)
{
  savvy::site_info variant;
  savvy::compressed_vector<float> genotypes;
  while (out && in.read(variant, genotypes))
  {
    if (update_info)
      savvy::update_info_fields(variant, genotypes, data_format);
    out.write(variant, genotypes);
  }

  if (regions.size())
  {
    for (auto it = regions.begin() + 1; it != regions.end(); ++it)
    {
      in.reset_bounds(*it);
      while (out && in.read(variant, genotypes))
      {
        if (update_info)
          savvy::update_info_fields(variant, genotypes, data_format);
        out.write(variant, genotypes);
      }
    }
  }

  if (out.fail())
    std::cerr << "Write Failure: Does input file have mixed ploidy?" << std::endl;

  return out.good() && !in.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}

int import_records(savvy::vcf::reader<1>& in, const std::vector<savvy::genomic_region>& regions, savvy::fmt data_format, bool update_info, savvy::sav::writer& out)
{
  // TODO: support regions without index.
  savvy::site_info variant;
  savvy::compressed_vector<float> genotypes;
  while (out && in.read(variant, genotypes))
  {
    if (update_info)
      savvy::update_info_fields(variant, genotypes, data_format);
    out.write(variant, genotypes);
  }

  if (out.fail())
    std::cerr << "Write Failure: Does input file have mixed ploidy?" << std::endl;
  return out.good() && !in.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <typename T>
int prep_reader_for_import(T& input, const import_prog_args& args)
{
  if (!input)
  {
    std::cerr << "Could not open file (" << args.input_path() << ")\n";
    return EXIT_FAILURE;
  }

  input.set_policy(args.empty_vector_policy());

  std::vector<std::string> sample_ids(input.samples().size());
  std::copy(input.samples().begin(), input.samples().end(), sample_ids.begin());
  if (args.subset_ids().size())
    sample_ids = input.subset_samples({args.subset_ids().begin(), args.subset_ids().end()});

  if (input.good())
  {
    auto headers = input.headers();
    headers.reserve(headers.size() + 3);
    headers.insert(headers.begin(), {"INFO","<ID=FILTER,Description=\"Variant filter\">"});
    headers.insert(headers.begin(), {"INFO","<ID=QUAL,Description=\"Variant quality\">"});
    headers.insert(headers.begin(), {"INFO","<ID=ID,Description=\"Variant ID\">"});

    savvy::sav::writer::options opts;
    opts.compression_level = args.compression_level();
    opts.block_size = args.block_size();
    if (args.index_path().size())
      opts.index_path = args.index_path();

    savvy::sav::writer output(args.output_path(), opts, sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end(), args.format());

    if (output.good())
    {
      if (args.sort_type())
      {
        //return (sort_and_write_records<std::vector<float>, less_than_comparator>((*args.sort_type()), input, args.format(), args.regions(), output, args.format(), args.update_info()) && !input.bad() ? EXIT_SUCCESS : EXIT_FAILURE);
      }
      else
      {
        return import_records(input, args.regions(), args.format(), args.update_info(), output);
      }
    }
  }

  return EXIT_FAILURE;
}

int import_main(int argc, char** argv)
{
  import_prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }


  if (args.regions().size())
  {
    savvy::vcf::indexed_reader<1> input(args.input_path(), args.regions().front(), args.bounding_point(), args.format());
    return prep_reader_for_import(input, args);
  }
  else
  {
    savvy::vcf::reader<1> input(args.input_path(), args.format());
    return prep_reader_for_import(input, args);
  }
}

int import_main2(int argc, char** argv)
{
  import_prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }


  if (args.regions().size())
  {
    fprintf(stderr, "Error: indexed import not yet implemented\n");
    return EXIT_FAILURE;
//    savvy::vcf::indexed_reader<1> input(args.input_path(), args.regions().front(), args.bounding_point(), args.format());
//    return prep_reader_for_import(input, args);
  }
  else
  {
    savvy::variant var;
    savvy::typed_value tmp_val;
    savvy::reader input(args.input_path());

    if (!input)
    {
      fprintf(stderr, "Error: could not open file %s\n", args.input_path().c_str());
      return EXIT_FAILURE;
    }

    auto hdrs = input.headers();

    auto gt_present = std::find_if(input.format_headers().begin(), input.format_headers().end(),
      [](const savvy::header_value_details& h) { return h.id == "GT"; }) != input.format_headers().end();

    bool remove_ph = input.phasing_status() != savvy::phasing::partial;
    if (gt_present && input.phasing_status() == savvy::phasing::unknown)
    {
      if (args.phasing() == savvy::phasing::unknown)
      {
        fprintf(stderr, "Error: phasing header not present, so --phasing must be specified\n");
        return EXIT_FAILURE;
      }

      std::string status;
      if (args.phasing() == savvy::phasing::none) status = "none";
      else if (args.phasing() == savvy::phasing::partial) remove_ph = false, status = "partial";
      else status = "full";

      hdrs.emplace_back("phasing", status);
    }

    if (args.pbwt_fields().size()) // TODO: check if PBWT headers already exist
    {
      hdrs.emplace_back("INFO", "<ID=_PBWT_RESET, Type=Flag>");
      for (const auto& f : args.pbwt_fields())
        hdrs.emplace_back("INFO", "<ID=_PBWT_SORT_" + f + ", Type=Flag, Format=" + f + ">");
    }

    if (input.phasing_status() == savvy::phasing::partial || args.phasing() == savvy::phasing::partial)
    {
      hdrs.emplace_back("FORMAT","<ID=PH, Type=Integer>");
    }

    // TODO: make args.subset_ids() a pointer to  allow for subsetting out all samples.
    auto subset_fn = args.subset_ids().empty() ? nullptr : ::savvy::detail::make_unique<savvy::typed_value::dense_subset_functor>(args.subset_ids(), input.samples());

    savvy::writer output(args.output_path(), savvy::file::format::sav2, hdrs, subset_fn ? subset_fn->id_intersection() : input.samples(), args.compression_level());
    output.set_block_size(args.block_size());

    std::size_t cnt = 0;
    while (output && input >> var)
    {
      if (remove_ph)
        var.set_format("PH", {});

      for (auto it = var.format_fields().begin(); it != var.format_fields().end(); ++it)
      {
        if (subset_fn)
        {
          it->second.copy_as_dense(tmp_val, *subset_fn);
          var.set_format(it->first, std::move(tmp_val));
        }

        if (args.sparse_fields().find(it->first) != args.sparse_fields().end())
        {
          it->second.copy_as_sparse(tmp_val);
          if (tmp_val.size() && static_cast<double>(tmp_val.non_zero_size()) / tmp_val.size() <= args.sparse_threshold())
            var.set_format(it->first, std::move(tmp_val)); // typed_value move operator implementation allows for reuse of tmp_val;
        }
        else if (it->second.is_sparse())
        {
          it->second.copy_as_dense(tmp_val);
          var.set_format(it->first, std::move(tmp_val));
        }
      }

      output << var;
      ++cnt;
    }

    return output.good() && !input.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
  }
}
#endif