# Hand-grading the newspaper claims

## Why you are doing this

An LLM (Llama 3.3 on Groq) reads all 1,324 scraped newspaper sentences and labels
each one. Nobody should take that on faith. So two people hand-label a random 80
of those same sentences, and we measure how often the humans and the LLM agree.

**Your labels are not training data.** Nothing is trained on them. They are a
report card. If the agreement number is good, we trust the LLM's labels on the
other ~1,244 claims and move on. If it's bad, we rewrite the rubric prompt in
`grade_claims.py` and re-run the LLM.

The number we compute is **Cohen's kappa**. It measures agreement *after
subtracting out the agreement you'd get by luck alone*. If 80% of claims are
"worsen," two people guessing "worsen" every time agree 80% of the time and have
learned nothing — kappa correctly scores that 0. Kappa ≥ 0.6 is "substantial"
and is the publishable floor; ≥ 0.8 is near-perfect.

Without this, the methods section reads "we asked an AI and believed it." That is
the first thing a judge will attack.

## The rules

1. **Two graders, independently.** Vincent takes one copy, Jeremy takes the other.
2. **Do not talk to each other while grading.** Not about individual claims, not
   about your general approach. If you discuss it, you've measured your
   conversation, not the reliability of the coding scheme.
3. **Do not look at `claims_graded.csv`.** That file has the LLM's answers in it.
   Peeking destroys the whole point.
4. **Grade only what the sentence itself says.** Not what you know happened next.
   You know 1929 was a catastrophe; the sentence doesn't. Judge the sentence.
5. **Disagreements are data, not failures.** The sentences you two split on are
   the ones the LLM is also getting wrong. Don't try to converge.

Budget 2–3 hours each. The sentences are short; the borderline cases are what
take the time, and those are the point.

## Setup

Each of you makes your own copy of the blank file:

```
cd JeremysShit/handgrade_newspapers
cp handgrade_BLANK.csv handgrade_vincent.csv
cp handgrade_BLANK.csv handgrade_jeremy.csv
```

Open your own copy in Excel or Google Sheets. Fill in the five `human_*` columns.
Leave everything else alone — `claim_id` is what joins your work back to the data.

The 80 claims are 8 from each of the 10 episodes, shuffled, so you cannot tell
which crisis (or control window) a sentence came from without looking at `date`.
That is deliberate. Try not to look at the date until after you've decided.

## What to type in each column

These come from OCR'd newspaper scans, so expect garbled text: `tho` for `the`,
`cor porations` for `corporations`, random `<yid` noise. Read through it.

### `human_is_prediction` → `yes` or `no`

`yes` only if the sentence makes a **falsifiable claim about future economic
conditions** — business, prices, employment, markets, prosperity, recession.

Say `no` for:
- retrospectives ("the panic ruined us")
- descriptions of the present ("stocks fell today")
- advertisements
- non-economic content (court rulings, land inspectors, elections)
- OCR garbage you can't parse

**The hard case, and the one that matters most.** A sentence like
*"President Roosevelt Thinks It Unnecessary as Panic Is Nearly Over"* is
partly a statement about the present and partly a claim that things will
improve. There is no right answer here — there is only your answer, and
whether the other grader and the LLM reach the same one. Make a call, write
your reasoning in `human_notes`, and move on.

Getting this boundary wrong in a *systematic* way is dangerous. If retrospectives
get filed as predictions, a paper printed the day after the 1929 crash saying
"we are ruined" scores as a brilliant forecast. That is fake signal, and it is
the single biggest threat to this project's credibility.

If `human_is_prediction` is `no`, **leave the next three columns blank.**

### `human_topic` → one of

| value | use when the claim is about |
|---|---|
| `general_business` | overall business conditions, prosperity, recession, "the outlook" |
| `prices` | inflation, deflation, cost of living, commodity prices |
| `employment` | jobs, unemployment, layoffs, hiring |
| `markets` | stocks, bonds, corporate earnings, Wall Street |
| `other` | economic but none of the above |

If it spans two, pick the one the sentence is *most* about.

### `human_direction` → one of

What does the claim say **economic conditions** will do?

- `improve` — recovery, prosperity returning, business picking up
- `worsen` — depression, panic, slump, downturn coming
- `no_change` — explicitly says things stay flat
- `unclear` — it's a prediction, but you genuinely can't tell which way

Careful with `prices` claims: *rising prices* is not automatically `improve`.
Ask what the sentence implies for **conditions overall**. "Prices will soar"
in a 1920 deflation context reads as recovery; the same phrase in 1945 reads
as ruinous inflation. Use the sentence's own framing, and if it gives you
nothing, use `unclear`.

### `human_confidence` → `assertive` or `hedged`

- `assertive` — *will*, *is certain*, *undoubtedly*, *there can be no doubt*
- `hedged` — *may*, *might*, *likely*, *is expected*, *if*, *some believe*

Judge the words, not the speaker's authority.

### `human_notes` → free text, optional

Use it whenever you hesitated. Write *why*. These notes are what you'll use to
rewrite the rubric if kappa comes back low, and they become the "coding
decisions" paragraph of the methods section. Notes on the cases you found hard
are worth more than notes on the easy ones.

## When you're both done

From `JeremysShit/`:

```
python handgrade_newspapers/kappa.py \
    --graders handgrade_newspapers/handgrade_vincent.csv \
              handgrade_newspapers/handgrade_jeremy.csv
```

`--graded` defaults to `claims_graded.csv` in the arm root; pass it explicitly
only if the LLM output lives somewhere else.

You get two blocks of numbers.

**HUMAN vs HUMAN** — is the rubric coherent? If you two can't agree with each
other, the instructions are ambiguous and no model can fix that. Rewrite this
file, regrade, try again.

**HUMAN vs LLM** — can we trust the LLM on the claims nobody checked? This is
the number that goes in the paper.

If `direction` comes back under 0.6, **stop.** Don't run `score_claims.py`, don't
retrain anything. Read the `human_notes` on the claims you disagreed about, find
the ambiguity, sharpen `RUBRIC_PROMPT` in `grade_claims.py`, re-run the LLM, and
recompute. Building a scoring rubric on labels the coders themselves cannot
reproduce just launders the disagreement into a number with a decimal point.

Report the final kappas in the paper, per field, with n. Report them even if
they're mediocre — a stated κ = 0.58 is a limitation, an unstated one is a hole.
