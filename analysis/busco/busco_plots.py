        for s in skipped:
            print(f"    skipped: {s}")

    if "synteny" not in skip:
        print(">>> Synteny plots")
        if args.synteny_pairs:
            pairs = [tuple(x.split(":", 1)) for x in args.synteny_pairs.split(",") if ":" in x]
        else:
            pairs = default_pairs(assemblies)
        if args.synteny_haps:
            pairs = pairs + hap_pairs(assemblies)
        print(f"    {len(pairs)} pair(s)")
        made, skipped = plot_synteny(pairs, assemblies, fulltables, args.final_dir,
                                     os.path.join(args.out_dir, "synteny"),
                                     args.synteny_orientation, args.chrs_limit)
        for m in made:
            print("    " + m)
        for s in skipped:
            print(f"    skipped: {s}")

    print(f">>> Done. Output under {args.out_dir}")


if __name__ == "__main__":
    main()
