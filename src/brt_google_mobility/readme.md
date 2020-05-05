If you hit an error running this task, it is likely that today's data are not yet available.
Run the command:
`orderly_run("rtm_incoming_ecdc", parameters = list(date_offset=1))` instead.

The `date_offset` lets you specify how many days back from "today" you want to pull the data.