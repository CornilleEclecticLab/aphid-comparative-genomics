🧬 General Note on Genome Version Consistency for Gene Annotations

When working with gene annotations—especially for gene families such as olfactory (OR) and gustatory receptor (GR) genes—it is essential to ensure compatibility between the genome assembly version and the corresponding annotation files.

❗Why This Matters:

Inconsistencies between genome versions and annotations can lead to:
	•	Internal stop codons or truncated CDS sequences
	•	Scaffold name mismatches
	•	Apparent frame shifts or stretches of ambiguous bases (Ns)
	•	Misinterpretation of gene integrity

🔄 General Recommendation:

If gene annotations were generated on an earlier version of a genome and are to be used with a newer assembly:
	•	Do not manually extract genes without checking genome version alignment
	•	Use tools like Liftoff or Lifton to transfer annotations accurately between versions
	•	Retain traceability by documenting the origin of annotations (version, tool, and source)
	•	Check gene integrity after transfer (CDS length, presence of stop codons, etc.)

🧰 Suggested Tools:
	•	Lifton: An improved version of Liftoff for robust annotation liftover
	•	Liftoff: Widely used for annotation projection across assemblies
	•	Custom scripts or pipelines to validate integrity after transfer
